
local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/protodunevd/simparams.jsonnet';
local fcl_params = {
    use_dnnroi: true,
};

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/protodunevd/sim.jsonnet';
local sim = sim_maker(params, tools);

// wire pitch dir (bottom drift):
// plane: 0 pitch dir: (0 -0.866025 0.5)
// plane: 1 pitch dir: (0 0.866025 0.5)
// plane: 2 pitch dir: (0 0 1)

// wire pitch dir (top drift):
// plane: 0 pitch dir: (0 0.866026 0.5)
// plane: 1 pitch dir: (0 -0.866026 0.5)
// plane: 2 pitch dir: (0 0 1)

local thetaXZ = 60*wc.deg;
local trackLen = 100;

// local stubby_top = {
//   tail: wc.point(100, 100, 100, wc.cm),
//   head: wc.point(100*(1 + std.tan(thetaXZ)), 100, 100*(1+1), wc.cm),
// };


local stubby_top1 = {
  tail: wc.point(10, -300, 10, wc.cm),
  head: wc.point(10 + trackLen* std.tan(thetaXZ), -300, 10 + trackLen, wc.cm),
};

local stubby_top2 = {
  tail: wc.point(10, -200, 10, wc.cm),
  head: wc.point(10 + trackLen* std.tan(thetaXZ), -200, 10 + trackLen, wc.cm),
};

local stubby_top3 = {
  tail: wc.point(10, -100, 10, wc.cm),
  head: wc.point(10 + trackLen* std.tan(thetaXZ), -100, 10 + trackLen, wc.cm),
};


// local stubby_bottom = {
//   tail: wc.point(-100, 100, 100, wc.cm),
//   head: wc.point(-100*(1 + std.tan(thetaXZ)), 100, 100*(1+1), wc.cm),
// };

local stubby_bottom1 = {
  tail: wc.point(-10, -300, 10, wc.cm),
  head: wc.point(-10 - trackLen* std.tan(thetaXZ), -300, 10 + trackLen, wc.cm),
};

local stubby_bottom2 = {
tail: wc.point(-10, -200, 10, wc.cm),
head: wc.point(-10 - trackLen* std.tan(thetaXZ), -200, 10 + trackLen, wc.cm),
};

local stubby_bottom3 = {
  tail: wc.point(-10, -100, 10, wc.cm),
  head: wc.point(-10 - trackLen* std.tan(thetaXZ), -100, 10 + trackLen, wc.cm),
};

local tracklist = [

  // {
  //   time: 0 * wc.us,
  //   charge: -500, // 5000 e/mm
  //   ray: stubby_top,
  // },

  // {
  //   time: 0 * wc.us,
  //   charge: -500,
  //   ray: stubby_bottom,
  //   // ray: params.det.bounds,
  // },

  {
    time: 0 * wc.us,
    charge: -500,
    ray: stubby_top1,
  },
  {
    time: 0 * wc.us,
    charge: -500,
    ray: stubby_top2,
  },
  {
    time: 0 * wc.us,
    charge: -500,
    ray: stubby_top3,
  },

  {
    time: 0 * wc.us,
    charge: -500,
    ray: stubby_bottom1,
  },
  {
    time: 0 * wc.us,
    charge: -500,
    ray: stubby_bottom2,
  },
  {
    time: 0 * wc.us,
    charge: -500,
    ray: stubby_bottom3,
  },

];

// local depos = sim.tracks(tracklist, step=1.0 * wc.mm);
local depos = sim.tracks(tracklist, step=0.1 * wc.mm);

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes-1);
local anode_idents = [anode.data.ident for anode in tools.anodes];

// local output = 'wct-sim-ideal-sig.npz';
// local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
local bagger = sim.make_bagger();
// signal plus noise pipelines
local sn_pipes = sim.splusn_pipelines;
// local analog_pipes = sim.analog_pipelines;

local perfect = import 'pgrapher/experiment/protodunevd/chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];

local nf_maker = import 'pgrapher/experiment/protodunevd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_override = if fcl_params.use_dnnroi then
{
    sparse: true,
    use_roi_debug_mode: true,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2],
    mp_tick_resolution: 4,
} else {
    sparse: true,
    use_multi_plane_protection: true,
    mp_tick_resolution: 4,
};

local sp_maker = import 'pgrapher/experiment/protodunevd/sp.jsonnet';
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local ts = {
    type: "TorchService",
    name: "dnnroi",
    data: {
        // model: "ts-model/unet-l23-cosmic500-e50.ts",
        model: "ts-model/CP49.ts",
        device: "cpu", // "gpucpu",
        concurrency: 1,
    },
};


local magoutput = 'protodunevd-sim-check-dnnsp.root';
local magnify = import 'pgrapher/experiment/protodunevd/magnify-sinks.jsonnet';
local magnifyio = magnify(tools, magoutput);

// Note: better switch to layers
local dnnroi = import 'pgrapher/experiment/protodunevd/dnnroi.jsonnet';

// local dnn_retaggers = [ g.pnode({
//     type: 'Retagger',
//     name: 'dnnsp-retag%d' %n,
//     data: {
//         // Note: retagger keeps tag_rules in an array to be like frame
//         // fanin/fanout but only uses first element.
//         tag_rules: [{
//             frame: {
//                 '.*': 'dnnsp',
//             },
//             merge: {
//                 'dnnsp*': 'dnnsp%d' %n,
//             },
//         }],
//     },
// }, nin=1, nout=1)
// for n in std.range(0, std.length(tools.anodes) - 1)
// ];

local parallel_pipes = [
  g.pipeline([
               sn_pipes[n],
               // magnifyio.orig_pipe[n],
               nf_pipes[n],
               magnifyio.raw_pipe[n],
               sp_pipes[n],
               // magnifyio.decon_pipe[n],
               magnifyio.debug_pipe[n],
             ] + if fcl_params.use_dnnroi then [
               dnnroi(tools.anodes[n], ts, output_scale=1.2),
               // dnn_retaggers[n],
               magnifyio.dnndecon_pipe[n],
             ] else [],
             'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];
local outtags = ['dnnsp%d' % n for n in std.range(0, std.length(tools.anodes) - 1)];
local parallel_graph = f.fanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf', outtags);



// // fixme: see https://github.com/WireCell/wire-cell-gen/issues/29
// local mega_anode = {
//   type: 'MegaAnodePlane',
//   name: 'meganodes',
//   data: {
//     anodes_tn: [wc.tn(anode) for anode in tools.anodes],
//   },
// };
// local make_noise_model = function(anode, csdb=null) {
//     type: "EmpiricalNoiseModel",
//     name: "empericalnoise-" + anode.name,
//     data: {
//         anode: wc.tn(anode),
//         // dft: wc.tn(tools.dft),
//         chanstat: if std.type(csdb) == "null" then "" else wc.tn(csdb),
//         spectra_file: params.files.noise,
//         nsamples: params.daq.nticks,
//         period: params.daq.tick,
//         wire_length_scale: 1.0*wc.cm, // optimization binning
//     },
//     // uses: [anode, tools.dft] + if std.type(csdb) == "null" then [] else [csdb],
//     uses: [anode] + if std.type(csdb) == "null" then [] else [csdb],
// };
// local noise_model = make_noise_model(mega_anode);
// local add_noise = function(model, n) g.pnode({
//     type: "AddNoise",
//     name: "addnoise%d-" %n + model.name,
//     data: {
//         rng: wc.tn(tools.random),
//         // dft: wc.tn(tools.dft),
//         model: wc.tn(model),
//   nsamples: params.daq.nticks,
//         replacement_percentage: 0.02, // random optimization
//     // }}, nin=1, nout=1, uses=[tools.random, tools.dft, model]);
//     }}, nin=1, nout=1, uses=[tools.random, model]);
// local noises = [add_noise(noise_model, n) for n in std.range(0,7)];
// 
// // local digitizer = sim.digitizer(mega_anode, name="digitizer", tag="orig");
// // "AnodePlane:anode110"
// // "AnodePlane:anode120"
// // "AnodePlane:anode111"
// // "AnodePlane:anode121"
// local digitizers = [
//     sim.digitizer(mega_anode, name="digitizer%d-" %n + mega_anode.name, tag="orig%d"%n)
//     for n in std.range(0,7)];
// 
// local frame_summers = [
//     g.pnode({
//         type: 'FrameSummer',
//         name: 'framesummer%d' %n,
//         data: {
//             align: true,
//             offset: 0.0*wc.s,
//         },
//     }, nin=2, nout=1) for n in std.range(0, 7)];
// 
// local actpipes = [g.pipeline([noises[n], digitizers[n]], name="noise-digitizer%d" %n) for n in std.range(0,7)];
// local util = import 'pgrapher/experiment/protodunevd/funcs.jsonnet';
// local pipe_reducer = util.fansummer('DepoSetFanout', analog_pipes, frame_summers, actpipes, 'FrameFanin');
// 
// 
// local magoutput = 'protodunevd-sim-check-dnnsp.root';
// local magnify = import 'pgrapher/experiment/protodunevd/magnify-sinks.jsonnet';
// local magnifyio = magnify(tools, magoutput);
// 
// // local fansel = g.pnode({
// //     type: "ChannelSplitter",
// //     name: "peranode",
// //     data: {
// //         anodes: [wc.tn(a) for a in tools.anodes],
// //         tag_rules: [{
// //             frame: {
// //                 '.*': 'orig%d'%ind,
// //             },
// //         } for ind in anode_idents/*anode_iota*/],
// //     },
// // }, nin=1, nout=nanodes, uses=tools.anodes);
// 
// local chsel = [
//   g.pnode({
//     type: 'ChannelSelector',
//     name: 'chsel%d' % n,
//     data: {
//       channels: util.anode_channels(n),
//       // tags: ['orig%d' % n], // traces tag
//     },
//   }, nin=1, nout=1)
//   for n in std.range(0, std.length(tools.anodes) - 1)
// ];
// 
// // Note: better switch to layers
// local dnnroi = import 'pgrapher/experiment/protodunevd/dnnroi.jsonnet';
// 
// local pipelines = [
//     g.pipeline([
//         chsel[n],
//         // magnifyio.orig_pipe[n],
// 
//         // nf_pipes[n],
//         // magnifyio.raw_pipe[n],
// 
//         sp_pipes[n],
//         // magnifyio.decon_pipe[n],
//         // magnifyio.threshold_pipe[n],
//         // magnifyio.debug_pipe[n], // use_roi_debug_mode=true in sp.jsonnet
//     ] + if fcl_params.use_dnnroi then [
//       dnnroi(tools.anodes[n], ts, output_scale=1.2),
//       magnifyio.dnndecon_pipe[n],
//     ] else [],
//                'nfsp_pipe_%d' % n)
//     for n in anode_iota
//     ];
// 
// // local fanpipe = f.fanpipe('FrameFanout', pipelines, 'FrameFanin', 'sn_mag_nf');
// local fanout_tag_rules = [
//           {
//             frame: {
//               '.*': 'orig%d' % tools.anodes[n].data.ident,
//             },
//             trace: {
//               // fake doing Nmult SP pipelines
//               //orig: ['wiener', 'gauss'],
//               //'.*': 'orig',
//             },
//           }
//           for n in std.range(0, std.length(tools.anodes) - 1)
//         ];
// 
// local anode_ident = [tools.anodes[n].data.ident for n in std.range(0, std.length(tools.anodes) - 1)];
// local fanin_tag_rules = [
//           {
//             frame: {
//               //['number%d' % n]: ['output%d' % n, 'output'],
//               '.*': 'framefanin',
//             },
//             trace: {
//               ['gauss%d'%ind]:'gauss%d'%ind,
//               ['wiener%d'%ind]:'wiener%d'%ind,
//               ['threshold%d'%ind]:'threshold%d'%ind,
//               ['dnnsp%d'%ind]:'dnnsp%d'%ind,
//             },
// 
//           }
//           for ind in anode_ident
//         ];
// local fanpipe = util.fanpipe('FrameFanout', pipelines, 'FrameFanin', 'nfsp', [], fanout_tag_rules, fanin_tag_rules);
// 
// // local fanin = g.pnode({
// //     type: 'FrameFanin',
// //     name: 'sigmerge',
// //     data: {
// //         multiplicity: nanodes,
// //         tags: [],
// //         tag_rules: [{
// //             trace: {
// //                 ['gauss%d'%ind]:'gauss%d'%ind,
// //                 ['wiener%d'%ind]:'wiener%d'%ind,
// //                 ['threshold%d'%ind]:'threshold%d'%ind,
// //             },
// //         } for ind in anode_iota],
// //     },
// // }, nin=nanodes, nout=1);
// 
// local retagger = g.pnode({
//     type: 'Retagger',
//     data: {
//         // Note: retagger keeps tag_rules in an array to be like frame
//         // fanin/fanout but only uses first element.
//         tag_rules: [{
//             // Retagger also handles "frame" and "trace" like
//             // fanin/fanout and also "merge" separately all traces
//             // like gaussN to gauss.
//             frame: {
//                 '.*': 'retagger',
//             },
//             merge: {
//                 'gauss\\d': 'gauss',
//                 'wiener\\d': 'wiener',
//                 'dnnsp\\d': 'dnnsp',
//             },
//         }],
//     },
// }, nin=1, nout=1);
// 
// // local fanpipe = g.intern(innodes=[fansel],
// //                          outnodes=[fanin],
// //                          centernodes=pipelines,
// //                          edges=
// //                          [g.edge(fansel, pipelines[n], n, 0) for n in anode_iota] +
// //                          [g.edge(pipelines[n], fanin, 0, n) for n in anode_iota],
// //                          name="fanpipe");

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;

// local graph = g.pipeline([depos, drifter, bagger, pipe_reducer, fanpipe, retagger, sink]);
local graph = g.pipeline([depos, drifter, bagger, parallel_graph, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot", "WireCellPytorch"],
        apps: ["Pgrapher"]
    }
};


// Finally, the configuration sequence which is emitted.

[cmdline] + g.uses(graph) + [app]
