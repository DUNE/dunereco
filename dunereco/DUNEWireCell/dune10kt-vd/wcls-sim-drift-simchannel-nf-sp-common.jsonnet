// Shared building blocks for the dune10kt-vd WC/LS
// sim-drift-simchannel-nf-sp job.  This "common" layer is detector-level glue:
// given the already-built `params`, `tools` and a resolved `fcl_params` object,
// it returns the reusable pnodes and node factories that BOTH the "single" and
// "multi" modal layers assemble into a graph.
//
// It does NOT resolve std.extVar() (that is the high-layer's job, see
// wcls-sim-drift-simchannel-nf-sp.jsonnet) and it does NOT build the final
// graph (that is the modal layers' job, see the -single / -multi files).
//
// Contract: the caller passes
//   params      : the fully-resolved dune10kt-vd parameter object
//   tools       : tools_maker(params) ALREADY narrowed to the anodes to process
//                 (one anode for single-*, a list for multi-*, all for full)
//   fcl_params  : object of the scalar fcl knobs (process_mode, use_dnnroi, ...)

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

local sim_maker = import 'pgrapher/experiment/dune10kt-vd/sim.jsonnet';
local nf_maker = import 'pgrapher/experiment/dune10kt-vd/nf.jsonnet';
local sp_maker = import 'pgrapher/experiment/dune10kt-vd/sp.jsonnet';
local chndb_perfect = import 'pgrapher/experiment/dune10kt-1x2x6/chndb-perfect.jsonnet';
local dnnroi = import 'pgrapher/experiment/dune-vd/dnnroi.jsonnet';

function(params, tools, fcl_params)

// WirePlaneLayer_t -> geo::_plane_proj.  Lexical local so the nested saver
// `data` objects below can reference it (inside those, `self` would bind to the
// data object, not to this common object).
local planemap_ = {
  dunevd_3view: { '1': 0, '2': 3, '4': 2 },
  default: { '1': 0, '2': 1, '4': 2 },
}[fcl_params.geo_planeid_labels];

{

  local sim = sim_maker(params, tools),

  // positional index over the (already narrowed) anode list
  anode_iota:: std.range(0, std.length(tools.anodes) - 1),

  // The input: an art-backed depo-set source feeding the drift stage.
  wcls_input:: {
    deposet: g.pnode({
      type: 'wclsSimDepoSetSource',
      name: '',
      data: {
        model: '',
        scale: -1,  // corrects a sign error in the SimDepoSource converter
        art_tag: fcl_params.input_depos_tag,  // upstream art producer of depos
        assn_art_tag: '',
        id_is_track: false,
      },
    }, nin=0, nout=1),
  },

  // WirePlaneLayer_t -> geo::_plane_proj (also exposed as a field).
  planemap:: planemap_,

  // A MegaAnodePlane spanning an arbitrary set of anodes.  The "single"/"full"
  // modes use one over all processed anodes; the "multi" mode builds one per
  // TPC so each per-TPC FrameSaver maps only its own channels.
  mega_anode_for(anodes):: {
    type: 'MegaAnodePlane',
    name: 'meganodes_' + std.join('_', [std.toString(a.data.ident) for a in anodes]),
    data: { anodes_tn: [wc.tn(a) for a in anodes] },
  },
  mega_anode:: self.mega_anode_for(tools.anodes),

  // Drift: a single DepoSetDrifter wrapping the per-depo Drifter.
  local drifter = sim.drifter,
  setdrifter:: g.pnode({
    type: 'DepoSetDrifter',
    data: { drifter: 'Drifter' },
  }, nin=1, nout=1, uses=[drifter]),

  // Per-anode simulation (signal+noise) pipelines (ductor->reframer->noise->
  // digitizer).  Each ends in a frame tagged "orig<ident>".
  sn_pipes:: sim.splusn_pipelines,

  // Per-anode "perfect" channel-noise DB + noise-filter pipelines.
  local chndb = [{
    type: 'OmniChannelNoiseDB',
    name: 'ocndbperfect%d' % n,
    data: chndb_perfect(params, tools.anodes[n], tools.field, n),
    uses: [tools.anodes[n], tools.field],
  } for n in self.anode_iota],
  nf_pipes:: [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n)
              for n in self.anode_iota],

  // Per-anode signal processing (OmnibusSigProc).  Tags traces "gauss<ident>",
  // "wiener<ident>" (and dnn ROI debug tags).
  local sp_override = if fcl_params.use_dnnroi then {
    sparse: true,
    use_roi_debug_mode: true,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2],
  } else {
    sparse: true,
  },
  local sp = sp_maker(params, tools, sp_override),
  sp_pipes:: [sp.make_sigproc(a) for a in tools.anodes],

  // DNN-ROI torch service and per-anode applier.
  ts:: {
    type: 'TorchService',
    name: 'dnnroi',
    data: {
      model: 'ts-model/unet-l23-cosmic500-e50.ts',
      device: 'gpucpu',
      concurrency: 1,
    },
  },
  dnnroi_node(n):: dnnroi(tools.anodes[n], self.ts, output_scale=1.2),

  // The wcls depo-flux writer (writes SimChannels to the art::Event and passes
  // the depo set through).  Three parameters are mode-dependent:
  //   simchan_label : the art product INSTANCE name of the std::vector<SimChannel>
  //   name          : the WCT node name == the art-output-converter name that
  //                   must appear in the module's `outputers` as
  //                   "wclsDepoFluxWriter:<name>"
  //   anodes        : the anode OBJECTS whose channels this writer tabulates.
  //                   The writer initializes one (dense) SimChannel per channel
  //                   of these anodes and projects only the depos flowing into
  //                   this node that fall inside their sensitive volumes.  Pass a
  //                   single-anode list to partition SimChannels per TPC (one
  //                   product per TPC); pass all processed anodes for one merged
  //                   product.  [default: all processed anodes, node "postdrift"]
  // `sparse`: false (default) pre-initializes one dense SimChannel per channel of
  // `anodes` EVERY event (visit() does this regardless of whether the node ran),
  // so the product is always nchan-sized -- even for a TPC the hydra busted.  Set
  // sparse=true to emit SimChannels ONLY for channels that actually received
  // charge (LArSoft-conventional); then a busted (no-depo) TPC yields a 0-size
  // product, matching the RawDigit/recob::Wire sparseness.
  depoflux_writer(simchan_label, name='postdrift', anodes=tools.anodes, sparse=false):: g.pnode({
    type: 'wclsDepoFluxWriter',
    name: name,
    data: {
      anodes: [wc.tn(anode) for anode in anodes],
      field_response: wc.tn(tools.field),
      tick: params.daq.tick,
      window_start: 0,
      window_duration: params.daq.readout_time,
      nsigma: 3.0,
      reference_time: 0,
      simchan_label: simchan_label,
      sed_label: fcl_params.input_depos_tag,
      sparse: sparse,
    },
  }, nin=1, nout=1, uses=anodes + [tools.field]),

  // ---- FrameSaver factories (the art output converters) -------------------
  // Each frame_tag in `tags` becomes an art product *instance* name; the art
  // module *label* is the FHiCL producer that hosts the WireCellToolkit module.
  // digitize=true -> raw::RawDigit, digitize=false -> recob::Wire.
  //
  // `anode` is the anode OBJECT (not tn) so it lands in the node's `uses`.

  // RawDigit saver (post-sim ADC).
  raw_saver(name, anode, tags):: g.pnode({
    type: 'wclsFrameSaver',
    name: name,
    data: {
      anode: wc.tn(anode),
      digitize: true,  // raw::RawDigit
      frame_tags: tags,
      pedestal_mean: 'native',
    },
  }, nin=1, nout=1, uses=[anode]),

  // recob::Wire saver (post-SP signals).  sparse=true (trick 2) writes the
  // zero-suppressed sigproc ROIs as sparse recob::Wire, so the product holds
  // only the non-zero regions instead of a dense all-ticks waveform -- far less
  // memory for the art::Event fan-in.
  sig_saver(name, anode, tags):: g.pnode({
    type: 'wclsFrameSaver',
    name: name,
    data: {
      plane_map: planemap_,
      anode: wc.tn(anode),
      digitize: false,  // recob::Wire
      sparse: true,     // zero-suppressed ROIs, not dense waveforms
      frame_tags: tags,
      frame_scale: [0.005 for t in tags],
      chanmaskmaps: [],
      nticks: params.daq.nticks,
    },
  }, nin=1, nout=1, uses=[anode]),

  // A terminal sink for frames whose content has already been saved to art.
  dumpframes(name):: g.pnode({
    type: 'DumpFrames',
    name: name,
    data: {},
  }, nin=1, nout=0),
}
