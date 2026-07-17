// This is a main entry point for configuring a wire-cell CLI job to
// simulate protoDUNE-SP.  It is simplest signal-only simulation with
// one set of nominal field response function.  It excludes noise.
// The kinematics are a mixture of Ar39 "blips" and some ideal,
// straight-line MIP tracks.
//
// Output is a Python numpy .npz file.

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/experiment/dune-vd/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';
local hs = import "pgrapher/common/helpers.jsonnet";

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params_maker = import 'pgrapher/experiment/dune-vd/params.jsonnet';
local response_plane = std.extVar('response_plane')*wc.cm;
local fcl_params = {
    G4RefTime: std.extVar('G4RefTime') * wc.us,
    response_plane: std.extVar('response_plane')*wc.cm,
    nticks: std.extVar('nticks'),
    ncrm: std.extVar('ncrm'),
    use_dnnroi: std.extVar('use_dnnroi'),
    process_crm: std.extVar('process_crm'),
};
local params = params_maker(fcl_params) {
  lar: super.lar {
    // Longitudinal diffusion constant
    DL: std.extVar('DL') * wc.cm2 / wc.ns,
    // Transverse diffusion constant
    DT: std.extVar('DT') * wc.cm2 / wc.ns,
    // Electron lifetime
    lifetime: std.extVar('lifetime') * wc.us,
    // Electron drift speed, assumes a certain applied E-field
    drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
  },
  files: super.files {
      wires: std.extVar('files_wires'),
      fields: [ std.extVar('files_fields'), ],
      noise: std.extVar('files_noise'),
  },
};

local tools_all = tools_maker(params);
local tools =
if fcl_params.process_crm == "partial"
then tools_all {anodes: [tools_all.anodes[n] for n in std.range(32, 79)]}
else if fcl_params.process_crm == "test1"
then tools_all {anodes: [tools_all.anodes[n] for n in [36]]}
else if fcl_params.process_crm == "test2"
then tools_all {anodes: [tools_all.anodes[n] for n in [36, 44]]}
else tools_all;

local sim_maker = import 'pgrapher/experiment/dune-vd/sim.jsonnet';
local sim = sim_maker(params, tools);

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);


local output = 'wct-sim-ideal-sig.npz';


//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
// local depos = sim.tracks(tracklist, step=1.0 * wc.mm);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);
local wcls_input = {
  depos: wcls.input.depos(name='', art_tag='IonAndScint'),
  // depos: wcls.input.depos(name='electron'),  // default art_tag="blopper"
  // DepoSet flavor input: needed by wclsDepoFluxWriter so IDepo::id() is the
  // SimEnergyDeposit index (id_is_track: false) for correct trackID mapping.
  deposet: g.pnode({
    type: 'wclsSimDepoSetSource',
    name: "",
    data: {
      model: "",
      scale: -1, //scale is -1 to correct a sign error in the SimDepoSource converter.
      art_tag: "IonAndScint", //name of upstream art producer of depos "label:instance:processName"
      id_is_track: false,
      assn_art_tag: "",
    },
  }, nin=0, nout=1),
};

// Collect all the wc/ls output converters for use below.  Note the
// "name" MUST match what is used in theh "outputers" parameter in the
// FHiCL that loads this file.
local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};

// WirePlaneLayer_t -> geo::_plane_proj
// U, V, W (1, 2, 4) -> U, V, W, Y (0, 1, 2, 3)
local planemaps = {
 dunevd_3view: {"1":0, "2":3, "4":2},
 default: {"1":0, "2":1, "4":2}
};
local planemap = planemaps[std.extVar("geo_planeid_labels")];

local wcls_output = {
  // ADC output from simulation
  // sim_digits: wcls.output.digits(name="simdigits", tags=["orig"]),
  sim_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'simdigits',
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['daq'],
      // nticks: params.daq.nticks,
      // chanmaskmaps: ['bad'],
      pedestal_mean: 'native',
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  // The output of signal processing.  Note, there are two signal
  // sets each created with its own filter.  The "gauss" one is best
  // for charge reconstruction, the "wiener" is best for S/N
  // separation.  Both are used in downstream WC code.
  // sp_signals: wcls.output.signals(name='spsignals', tags=['gauss', 'wiener']),
  sp_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'spsignals',
    data: {
      plane_map: planemap,
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['gauss', 'wiener','dnnsp'],
      frame_scale: [0.005, 0.005, 0.005],
      chanmaskmaps: [],
      nticks: params.daq.nticks,
    },
  }, nin=1, nout=1, uses=[mega_anode]),
};

//local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
local setdrifter = g.pnode({
            type: 'DepoSetDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1,
        uses=[drifter]);
local bagger = sim.make_bagger();
// local bagger = g.pnode({
//   type: 'DepoBagger',
//   name: 'bagger',
//   data: {
//     gate: [-250 * wc.us, 2750 * wc.us],  // fixed
//   },
// }, nin=1, nout=1);

// signal plus noise pipelines
// local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local perfect = import 'pgrapher/experiment/dune10kt-1x2x6/chndb-perfect.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n),
  uses: [tools.anodes[n], tools.field],  // pnode extension
} for n in anode_iota];

//local chndb_maker = import 'pgrapher/experiment/pdsp/chndb.jsonnet';
//local noise_epoch = "perfect";
//local noise_epoch = "after";
//local chndb_pipes = [chndb_maker(params, tools.anodes[n], tools.fields[n]).wct(noise_epoch)
//                for n in std.range(0, std.length(tools.anodes)-1)];
local nf_maker = import 'pgrapher/experiment/dune-vd/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/dune-vd/sp.jsonnet';
local sp_override = if fcl_params.use_dnnroi then
{
    sparse: true,
    use_roi_debug_mode: true,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2]
} else {
    sparse: true,
};
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

// ── DNN-ROI configuration (protodune style: explicit parameters) ────
// FD-VD uses the ProtoDUNE-VD 6-channel model (same CRP anode design).
// FP32 best KD (6-ch) PDVD model.  Also resolvable via WIRECELL_PATH as
// 'dnnroi/pdvd/pipe_distill_transformer_6ch.ts'; located in Hugging Face.
local dnnroi_model = '/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/WireCell/dune/dnn-roi/pdvd/20260615/pipe_distill_transformer_6ch.ts';
local dnnroi_device = 'cpu';                   // 'cpu' or 'gpu'
local dnnroi_nchan = 6;                        // 6 = production KD/QAT models
local dnnroi_concurrency = 1;
local dnnroi_nticks = fcl_params.nticks;
local dnnroi_tick_per_slice = 4;               // training rebin=4
local dnnroi_output_scale = 1.0;               // 6-ch .ts bakes normalization; no ad-hoc scale
local dnnroi_mask_thresh = 0.2;
local dnnroi_nchunks = 1;

local ts = {
    type: "TorchService",
    name: "dnnroi",
    data: {
        model: dnnroi_model,
        device: dnnroi_device,
        concurrency: dnnroi_concurrency,
    },
};

// Per-anode per-plane 6-channel DNN-ROI subgraph (protodune style).
local dnnroi = import 'pgrapher/experiment/dune-vd/dnnroi_pp.jsonnet';

// Truth labeling: wclsDepoFluxWriter (switched from wclsSimChannelSink).
local wcls_depoflux_writer = g.pnode({
  type: 'wclsDepoFluxWriter',
  name: 'postdrift',
  data: {
    anodes: [wc.tn(anode) for anode in tools.anodes],
    field_response: wc.tn(tools.field),
    tick: 0.5 * wc.us,
    window_start: 0.0 * wc.ms,
    window_duration: self.tick * fcl_params.nticks,
    nsigma: 3.0,
    reference_time: fcl_params.G4RefTime,
    energy: 0, // energy: 0 means use the true depo energy (equivalent to SimChannelSink use_energy=true)
    simchan_label: 'simpleSC',
    sed_label: 'IonAndScint',
    sparse: false,
    // Smear the truth SimChannel to reco (SP) resolution so the 2D truth labels
    // register onto the deconvolved reco charge (added in quadrature to the
    // post-drift depo diffusion). Values = dune-vd SP-filter recipe
    // (sp-filters.jsonnet):
    //   smear_long [ticks] = sigma_t/tick, sigma_t = 1/(2*pi*0.12MHz) = 1.326us
    //   smear_tran [pitch] = 1/(2*sqrt(pi)*k): Wire_ind k=0.75 (U,V); Wire_col k=3.0 (W)
    smear_long: 2.6526,
    smear_tran: [0.37612, 0.37612, 0.09403],
  },
}, nin=1, nout=1, uses=tools.anodes + [tools.field]);

local magoutput = 'mag-sim-sp.root';
local magnify = import 'pgrapher/experiment/dune-vd/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);

local multipass = [
  g.pipeline([
                // wcls_simchannel_sink[n],
                sn_pipes[n],
                // sinks.orig_pipe[n],
                // nf_pipes[n],
                sp_pipes[n],
                // sinks.decon_pipe[n],
                // sinks.debug_pipe[n], // use_roi_debug_mode=true in sp.jsonnet
             ] + if fcl_params.use_dnnroi then [
                 dnnroi(tools.anodes[n], ts,
                        nticks=dnnroi_nticks,
                        tick_per_slice=dnnroi_tick_per_slice,
                        output_scale=dnnroi_output_scale,
                        mask_thresh=dnnroi_mask_thresh,
                        nchunks=dnnroi_nchunks,
                        nchan=dnnroi_nchan),
                //  sinks.dnnroi_pipe[n],
             ] else [],
             'multipass%d' % n)
  for n in anode_iota
];

local f = import 'pgrapher/experiment/dune-vd/funcs.jsonnet';
// local outtags = ['gauss%d' % anode.data.ident for anode in tools.anodes];
local outtags = [];
local tag_rules = {
    frame: {
        '.*': 'framefanin',
    },
    trace: {['gauss%d' % anode.data.ident]: ['gauss%d' % anode.data.ident] for anode in tools.anodes}
        + {['wiener%d' % anode.data.ident]: ['wiener%d' % anode.data.ident] for anode in tools.anodes}
        + {['threshold%d' % anode.data.ident]: ['threshold%d' % anode.data.ident] for anode in tools.anodes}
        + {['dnnsp%d' % anode.data.ident]: ['dnnsp%d' % anode.data.ident] for anode in tools.anodes},
};
local bi_manifold =
    if fcl_params.ncrm == 36
    then f.multifanpipe('DepoSetFanout', multipass, 'FrameFanin', [1,6], [6,6], [1,6], [6,6], 'sn_mag', outtags, tag_rules)
    else if fcl_params.ncrm == 48 || fcl_params.process_crm == "partial"
    then f.multifanpipe('DepoSetFanout', multipass, 'FrameFanin', [1,8], [8,6], [1,8], [8,6], 'sn_mag', outtags, tag_rules)
    else if fcl_params.process_crm == "test1"
    then f.multifanpipe('DepoSetFanout', multipass, 'FrameFanin', [1,1], [1,1], [1,1], [1,1], 'sn_mag', outtags, tag_rules)
    else if fcl_params.process_crm == "test2"
    then f.multifanpipe('DepoSetFanout', multipass, 'FrameFanin', [1,2], [2,1], [1,2], [2,1], 'sn_mag', outtags, tag_rules)
    else if fcl_params.ncrm == 112
    then f.multifanpipe('DepoSetFanout', multipass, 'FrameFanin', [1,8,16], [8,2,7], [1,8,16], [8,2,7], 'sn_mag', outtags, tag_rules);

local retagger = g.pnode({
  type: 'Retagger',
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like gaussN to gauss.
      frame: {
        '.*': 'retagger',
      },
      merge: {
        'gauss\\d+': 'gauss',
        'wiener\\d+': 'wiener',
        'dnnsp\\d+': 'dnnsp',
      },
    }],
  },
}, nin=1, nout=1);

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;

local graph = g.pipeline([wcls_input.deposet, setdrifter, wcls_depoflux_writer, bi_manifold, retagger, wcls_output.sp_signals, sink]);

local app = {
    type: std.extVar('engine'), //Pgrapher, TbbFlow
    data: {
        edges: g.edges(graph),
    },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
