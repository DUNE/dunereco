// This is a main entry point for configuring a wire-cell CLI job to
// simulate protoDUNE-SP.  It is simplest signal-only simulation with
// one set of nominal field response function.  It excludes noise.
// The kinematics are a mixture of Ar39 "blips" and some ideal,
// straight-line MIP tracks.
//
// Output is a Python numpy .npz file.

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';
// Per-anode per-plane 6-channel DNN-ROI subgraph (protodune style).
local dnnroi = import 'pgrapher/experiment/dune10kt-vd/dnnroi_pp.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local response_plane = std.extVar('response_plane')*wc.cm;
local fcl_params = {
    G4RefTime: std.extVar('G4RefTime') * wc.us,
    response_plane: std.extVar('response_plane')*wc.cm,
    nticks: std.extVar('nticks'),
    ncrm: std.extVar('ncrm'),
    use_dnnroi: std.extVar('use_dnnroi'),
    process_crm: std.extVar('process_crm'),
    use_hydra: std.extVar('use_hydra'),
    save_rawdigits: std.extVar('save_rawdigits'),
    adc_resolution: std.extVar('adc_resolution'),
};
local params_maker =
if fcl_params.ncrm == 320 then import 'pgrapher/experiment/dune10kt-vd/params-10kt.jsonnet'
else import 'pgrapher/experiment/dune10kt-vd/params.jsonnet';
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
then tools_all {anodes: [tools_all.anodes[n] for n in [5]]}
else if fcl_params.process_crm == "test2"
then tools_all {anodes: [tools_all.anodes[n] for n in [0,1,4,5]]}
else tools_all;

local sim_maker = import 'pgrapher/experiment/dune10kt-vd/sim.jsonnet';
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
    deposet: g.pnode({
            type: 'wclsSimDepoSetSource',
            name: "", 
            data: {
                model: "", 
                scale: -1, //scale is -1 to correct a sign error in the SimDepoSource converter.
                art_tag: "IonAndScint", //name of upstream art producer of depos "label:instance:processName"
                id_is_track: false, // IDepo::id() is the SimEnergyDeposit index, needed by wclsDepoFluxWriter
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
local nf_maker = import 'pgrapher/experiment/dune10kt-vd/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/dune10kt-vd/sp.jsonnet';
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
    // post-drift depo diffusion). Values = dune10kt-vd SP-filter recipe
    // (sp-filters.jsonnet):
    //   smear_long [ticks] = sigma_t/tick, sigma_t = 1/(2*pi*0.12MHz) = 1.326us
    //   smear_tran [pitch] = 1/(2*sqrt(pi)*k): Wire_ind k=0.75 (U,V); Wire_col k=3.0 (W)
    smear_long: 2.6526,
    smear_tran: [0.37612, 0.37612, 0.09403],
  },
}, nin=1, nout=1, uses=tools.anodes + [tools.field]);

local magoutput = 'mag.root';
local magnify = import 'pgrapher/experiment/dune-vd/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);


local full_sim_pipes = [
  g.pipeline([
                sn_pipes[n],
                // sinks.orig_pipe[n],
             ],
             'multipass%d' % n)
  for n in anode_iota
];

local full_sp_pipes = [
  g.pipeline([
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
             'full_sp_pipes%d' % n)
  for n in anode_iota
];

local multipass = [
  g.pipeline([
                // wcls_simchannel_sink[n],
                sn_pipes[n],
                // sinks.orig_pipe[n],
                // nf_pipes[n],
                sp_pipes[n],
                sinks.decon_pipe[n],
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

local f = import 'pgrapher/common/funcs.jsonnet';
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


local make_switch_pipe = function(sim, sp, anode ) {
    local d2f = g.pipeline([sim, sp]),
    local ds_filter = g.pnode({
        type: "DepoSetFilter",
        name: "ds-filter-switch-%d" % anode.data.ident,
        data: {anode: wc.tn(anode)},
        }, nin=1, nout=1, uses=[anode]),
    local dorb = g.pnode({
        type: "DeposOrBust",
        name: "dorb-switch-%d" % anode.data.ident,
        }, nin=1, nout=2),
    local frame_sync = g.pnode({
        type: "FrameSync",
        name: "frame-sync-switch-%d" % anode.data.ident,
        }, nin=2, nout=1),
    // direct pipe
    local ret1 = g.pipeline([ds_filter, d2f]),
    // hydra shortcut
    local ret2 =  g.intern(
        innodes=[ds_filter],
        outnodes=[frame_sync],
        centernodes=[dorb, d2f],
        edges=
            [g.edge(ds_filter, dorb, 0, 0),
            g.edge(dorb, d2f, 0, 0),
            g.edge(d2f, frame_sync, 0, 0),
            g.edge(dorb, frame_sync, 1, 1)]),
    // special case to tapout and sync rawdigits
    local fout_bust = g.pnode({
        type: "FrameFanout",
        name: "fout-switch-bust-%d" % anode.data.ident,
        data: {multiplicity: 2},
        }, nin=1, nout=2),
    local fout_rawdigits = g.pnode({
        type: "FrameFanout",
        name: "fout-switch-rawdigits-%d" % anode.data.ident,
        data: {multiplicity: 2},
        }, nin=1, nout=2),
    local frame_sync_rawdigits = g.pnode({
        type: "FrameSync",
        name: "frame-sync-switch-rawdigits-%d" % anode.data.ident,
        }, nin=2, nout=1),
    
    local dump_rawdigits = g.pnode(
        { type: 'DumpFrames', name:"switch-dump-rawdigits-%d" % anode.data.ident, data:{} },
        nin=1, nout=0),
    // hydra shortcut with rawdigits
    local ret3 = g.intern(
        innodes=[ds_filter],
        outnodes=[frame_sync],
        centernodes=[dorb, sim, sp, fout_bust, fout_rawdigits, frame_sync_rawdigits, dump_rawdigits],
        edges=
            [g.edge(ds_filter, dorb, 0, 0),
            g.edge(dorb, sim, 0, 0),
            g.edge(dorb, fout_bust, 1, 0),
            g.edge(sim, fout_rawdigits, 0, 0),
            g.edge(fout_rawdigits, sp, 0, 0),
            g.edge(sp, frame_sync, 0, 0),
            g.edge(fout_bust, frame_sync, 0, 1),
            g.edge(fout_rawdigits, frame_sync_rawdigits, 1, 0),
            g.edge(fout_bust, frame_sync_rawdigits, 1, 1),
            g.edge(frame_sync_rawdigits, dump_rawdigits, 0, 0),
            ]
    ),
    ret: if fcl_params.save_rawdigits then ret3 else ret2,
}.ret;

local switch_pipes = [
    make_switch_pipe(full_sim_pipes[n], full_sp_pipes[n], tools.anodes[n]),
    for n in std.range(0, std.length(tools.anodes) - 1)
];

local process_pipes = if fcl_params.use_hydra then switch_pipes else multipass;

local bi_manifold =
    if fcl_params.process_crm == "test2"
    then f.multifanpipe('DepoSetFanout', process_pipes, 'FrameFanin', [1,4], [4,1], [1,4], [4,1], 'sn_mag', outtags, tag_rules)
    else f.multifanpipe('DepoSetFanout', process_pipes, 'FrameFanin', [1,2,8,32], [2,4,4,10], [1,2,8,32], [2,4,4,10], 'sn_mag', outtags, tag_rules);


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


// Build an incomplete subgraph ending to be spliced for saving out frames 
// local osimfanin = g.pnode({ 
//     type: 'FrameFanin',
//     name:"osimfanin",
//     data:{
//         multiplicity: std.length(tools.anodes),
//         tags: ['orig%d' % n for n in anode_iota],
//     } 
// }, nin=std.length(tools.anodes), nout=1);
local osimfanin = 
if fcl_params.process_crm == "test2"
then f.multifanin('FrameFanin', [1,4], [4,1], 'osimfanin', ['orig%d' % n for n in anode_iota])
else f.multifanin('FrameFanin', [1,2,8,32], [2,4,4,10], 'osimfanin', ['orig%d' % n for n in anode_iota]);
local osimdump = g.pnode({ type: 'DumpFrames', name:"osimdump", data:{} }, nin=1, nout=0);
local osimtagger = g.pnode({
  type: 'Retagger',
  name: 'osimtagger',
  data: {
    tag_rules: [{
      frame: {'.*': 'daq',},
      merge: {
        'orig\\d': 'orig',
      },
    }],
  },
}, nin=1, nout=1);
local osimsaver = g.pipeline([osimfanin, osimtagger, wcls_output.sim_digits, osimdump]);

local edge_selector(e) = std.startsWith(e.tail.node, "FrameSync:frame-sync-switch-rawdigits");
// local fanout_factory(n,e) = { type:'FrameFanout', name:"splice%d"%n };
local fanout_factory(n,e) = { type:'FrameFanout', name:"splice%d"%n, data:{multiplicity: 2} }; // "2-wire" splice



local main_graph = g.pipeline([wcls_input.deposet, setdrifter, wcls_depoflux_writer, bi_manifold, retagger, wcls_output.sp_signals, sink]);

local graph = g.splice(main_graph, osimsaver, edge_selector, fanout_factory);

local app = {
    type: 'TbbFlow', //Pgrapher, TbbFlow
    data: {
        edges: g.edges(graph),
    },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
