// WC/LS DNN-ROI entry point for ProtoDUNE-HD.
//
// This is the LArSoft (WC/LS) counterpart of the pure-WCT driver
// pgrapher/experiment/pdhd/wct-nf-sp-dnnroi.jsonnet.  It reads
// RawDigits from the art::Event (wclsRawFrameSource), runs NF + SP +
// DNN-ROI + L1SP-after-DNN per anode, fans the per-anode frames back
// together and saves the result via wclsFrameSaver:spsaver.
//
// Built from pdhd/wcls-nf-sp.jsonnet (the canonical WC/LS NF-SP entry)
// with the DNN-ROI + L1SP-after-DNN subgraph from wct-nf-sp-dnnroi.jsonnet
// spliced in per anode.
//
// Defaults (per deployment request):
//   - DNN-ROI enabled (use_dnnroi=true)
//   - inference on CPU (dnnroi_device / l1sp_pd_dnn_device = "cpu")
//   - L1SP-after-DNN enabled in 'hybrid' mode (l1sp_pd_mode='hybrid')
//
// The DNN-ROI input tags (loose_lf*, mp2_roi*, mp3_roi*, tight_lf*,
// decon_charge*, gauss*) require OmnibusSigProc to run in debug +
// multi-plane-protection mode, so sp_override forces those on whenever
// DNN-ROI is enabled.  L1SP stays OFF inside sp.make_sigproc (the SP
// auxiliary tags would be dropped through its FrameMerger); L1SP runs
// strictly AFTER DNN-ROI inside the l1sp_after_dnnroi envelope.

local epoch = std.extVar('epoch');  // eg "dynamic", "after", "before", "perfect"
local reality = std.extVar('reality');
local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"

// When true, ALSO save the traditional SP output (gauss/wiener) taken
// BEFORE DNN-ROI, in addition to the post-DNN result.  This is done by
// splicing a 2-way fanout onto every OmnibusSigProc output edge: one
// branch continues into the DNN-ROI / L1SP chain (unchanged), the other
// is fanned-in, retagged and written by wclsFrameSaver:tradspsaver.
// The matching 'tradspsaver' name must appear in the FHiCL 'outputers'.
local save_tradsp = true;


local wc = import 'wirecell.jsonnet';
local g = import 'pgraph.jsonnet';

local raw_input_label = std.extVar('raw_input_label');  // eg "daq"


local data_params = import 'params.jsonnet';
local simu_params = import 'simparams.jsonnet';
local params = if reality == 'data' then data_params else simu_params;


local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

// ── DNN-ROI / L1SP configuration ────────────────────────────────────
// DNN-ROI is on by default; pass use_dnnroi: "false" (string) in the
// FHiCL structs to fall back to bare NF+SP.
local use_dnnroi = std.extVar('use_dnnroi');

// FP32 best KD (6-ch) PDHD model.  Resolved via WIRECELL_PATH, located in Hugging face.
// local dnnroi_model = 'dnnroi/pdhd/pipe_distill_transformer_6ch.ts';
local dnnroi_model = '/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/WireCell/dune/dnn-roi/pdhd/20260615/pipe_distill_transformer_6ch.ts';
local dnnroi_device = 'cpu';                   // 'cpu' or 'gpu'.  INT8 graph is CPU-only.
local dnnroi_nchan = 6;                        // 6 = production KD/QAT models
local dnnroi_concurrency = 1;
local dnnroi_nticks = 6000;
local dnnroi_tick_per_slice = 4;               // training rebin=4
local dnnroi_output_scale = 1.0;
local dnnroi_mask_thresh = 0.2;
local dnnroi_nchunks = 1;

// SP ROI-refinement tunes (mirror wct-nf-sp-dnnroi.jsonnet defaults).
local apa0_w_roi_tune = true;
local sp_roi_mad_rms = true;
local sp_w_col_break_roi_tune = true;

// L1SP-after-DNN: enabled by default in 'hybrid' mode.
local use_l1sp_dnn = true;
local l1sp_pd_mode = 'hybrid';                 // 'process'|'dump'|'dnn'|'hybrid'|''
local l1sp_pd_adj_enable = true;
local l1sp_pd_adj_max_hops = 3;
local l1sp_pd_dump_path = '';
local l1sp_pd_wf_dump_path = '';
local l1sp_pd_dump_all_rois = false;
// L1SP DNN tagger (mode 'dnn'||'hybrid').  Model resolved via WIRECELL_PATH; located in Hugging face.
// local l1sp_pd_dnn_model = 'l1sp/pdhd/l1sp_dnn_pdhd_v1.ts';
local l1sp_pd_dnn_model = '/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/WireCell/dune/l1sp/pdhd/20260615/l1sp_dnn_pdhd_v1.ts';
local l1sp_pd_dnn_device = 'cpu';
local l1sp_pd_dnn_concurrency = 1;
local l1sp_pd_dnn_threshold = 0.5;
local l1sp_pd_dnn_window_ticks = 256;
local l1sp_pd_dnn_debug_path = '';
local l1sp_pd_adj_dnn_veto = true;
// Loose-heur pre-filter overrides (DNN-chain values).
local l1sp_pd_gmax_min = 1500.0;
local l1sp_pd_min_length = 30;
local l1sp_pd_energy_frac_thr = 0.66;

local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';
// DNN-ROI input tags require OmnibusSigProc debug + multi-plane-protection.
local sp_override = { sparse: sigoutform == 'sparse' }
                    + (if use_dnnroi
                       then { use_roi_debug_mode: true, use_multi_plane_protection: true }
                       else {});

local use_resampler = (reality == 'data');

// Collect the WC/LS input converters for use below.  Make sure the
// "name" argument matches what is used in the FHiCL that loads this
// file.  In particular if there is no ":" in the inputer then name
// must be the emtpy string.
local wcls_input = {
  adc_digits: g.pnode({
    type: 'wclsRawFrameSource',
    name: '',
    data: {
      art_tag: raw_input_label,
      frame_tags: ['orig'],  // this is a WCT designator
      tick: 512 * wc.ns,
      // nticks: params.daq.nticks,
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
local wcls_output = {
  // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'nfsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['raw'],
      chanmaskmaps: ['bad'],
    },
  }, nin=1, nout=1, uses=[mega_anode]),


  // The output of signal processing.  After DNN-ROI the "gauss" tag
  // carries the DNN-ROI output (U+V from the model, W passthrough),
  // optionally refined by L1SP; "wiener" carries the SP Wiener /
  // per-channel threshold summary.
  sp_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'spsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['gauss', 'wiener'],
      frame_scale: [0.001, 0.001],
      chanmaskmaps: [],
      nticks: -1,
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  // Traditional SP output (gauss/wiener) tapped BEFORE DNN-ROI.  Only
  // wired into the graph when save_tradsp=true (see the splice below).
  // Mirrors sp_signals but writes a distinct art product so the
  // pre-DNN and post-DNN SP can be compared event-by-event.
  tradsp_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'tradspsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      // Use distinct instance names so they don't collide with the
      // post-DNN spsaver's gauss/wiener art products.
      // No underscores — art rejects them in product instance names.
      frame_tags: ['tradspgauss', 'tradspwiener'],
      frame_scale: [0.001, 0.001],
      chanmaskmaps: [],
      nticks: -1,
    },
  }, nin=1, nout=1, uses=[mega_anode]),
};

local base = import 'chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: base(params, tools.anodes[n], tools.field, n) { dft: wc.tn(tools.dft) },
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in std.range(0, std.length(tools.anodes) - 1)];

local nf_maker = import 'pgrapher/experiment/pdhd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp = sp_maker(params, tools, sp_override);
// L1SP OFF (l1sp_pd_mode='') inside SP so the DNN-ROI debug input tags
// survive; L1SP runs after DNN-ROI in the envelope below.
local sp_pipes = [sp.make_sigproc(a,
                                  l1sp_pd_mode='',
                                  apa0_w_roi_tune=apa0_w_roi_tune,
                                  roi_mad_rms=sp_roi_mad_rms,
                                  w_col_break_roi_tune=sp_w_col_break_roi_tune)
                  for a in tools.anodes];

// TorchService shared by all per-anode DNN-ROI nodes.
local ts = {
  type: 'TorchService',
  name: 'dnnroi_pdhd',
  data: {
    model: dnnroi_model,
    device: dnnroi_device,
    concurrency: dnnroi_concurrency,
  },
};

local dnnroi_maker = import 'pgrapher/experiment/pdhd/dnnroi_pp.jsonnet';
local dnnroi_inner_pipes = [dnnroi_maker(tools.anodes[n], ts,
                                         nticks=dnnroi_nticks,
                                         tick_per_slice=dnnroi_tick_per_slice,
                                         output_scale=dnnroi_output_scale,
                                         mask_thresh=dnnroi_mask_thresh,
                                         nchunks=dnnroi_nchunks,
                                         nchan=dnnroi_nchan)
                            for n in std.range(0, std.length(tools.anodes) - 1)];

// L1SP-after-DNN envelope.
local l1sp_dnn_maker = import 'pgrapher/experiment/pdhd/l1sp_after_dnnroi.jsonnet';

// TorchService for the L1SP DNN tagger.  Built when mode is 'dnn' or
// 'hybrid'; null otherwise so the heuristic-only path doesn't pull
// libtorch a second time.
local l1sp_torch_service =
  if (l1sp_pd_mode == 'dnn' || l1sp_pd_mode == 'hybrid') then {
    type: 'TorchService',
    name: 'l1sp_dnn_pdhd',
    data: {
      model: l1sp_pd_dnn_model,
      device: l1sp_pd_dnn_device,
      concurrency: l1sp_pd_dnn_concurrency,
    },
  } else null;

// Rename the inner DNN-ROI subgraph's unified trace tag dnnsp%d -> gauss%d
// for the L1SP-off path so the downstream fan-in (which expects gauss%d)
// is satisfied.  Used only when use_l1sp_dnn=false.
local dnn_to_gauss_retag = function(n)
  g.pnode({
    type: 'Retagger',
    name: 'dnnsp_to_gauss%d' % n,
    data: {
      tag_rules: [{
        frame: { ['dnnsp%d' % n]: 'gauss%d' % n },
        merge: { ['dnnsp%d' % n]: 'gauss%d' % n },
      }],
    },
  }, nin=1, nout=1);

local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: std.range(2560 * n, 2560 * (n + 1) - 1),
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local resamplers_config = import 'pgrapher/common/resamplers.jsonnet';
local load_resamplers = resamplers_config(g, wc, tools);
local resamplers = load_resamplers.resamplers;

// Per-anode SP + DNN-ROI + (optional) L1SP-after-DNN segment.  Mirrors
// wct-nf-sp-dnnroi.jsonnet's sp_dnn_l1sp_segment.
local sp_dnn_l1sp_segment(n) =
  if use_dnnroi && use_l1sp_dnn
  then [l1sp_dnn_maker(tools.anodes[n], sp_pipes[n], dnnroi_inner_pipes[n],
                       tools, params,
                       l1sp_pd_adj_enable=l1sp_pd_adj_enable,
                       l1sp_pd_adj_max_hops=l1sp_pd_adj_max_hops,
                       l1sp_pd_dump_mode=l1sp_pd_mode,
                       l1sp_pd_dump_path=l1sp_pd_dump_path,
                       l1sp_pd_wf_dump_path=l1sp_pd_wf_dump_path,
                       l1sp_pd_dump_all_rois=l1sp_pd_dump_all_rois,
                       l1sp_pd_torch_service=l1sp_torch_service,
                       l1sp_pd_dnn_threshold=l1sp_pd_dnn_threshold,
                       l1sp_pd_adj_dnn_veto=l1sp_pd_adj_dnn_veto,
                       l1sp_pd_dnn_window_ticks=l1sp_pd_dnn_window_ticks,
                       l1sp_pd_dnn_debug_path=l1sp_pd_dnn_debug_path,
                       l1sp_pd_gmax_min=l1sp_pd_gmax_min,
                       l1sp_pd_min_length=l1sp_pd_min_length,
                       l1sp_pd_energy_frac_thr=l1sp_pd_energy_frac_thr)]
  else [sp_pipes[n]]
       + (if use_dnnroi
          then [dnnroi_inner_pipes[n], dnn_to_gauss_retag(tools.anodes[n].data.ident)]
          else []);

local nfsp_pipes = [
  g.pipeline(
    [chsel_pipes[n]]
    + (if use_resampler then [resamplers[n]] else [])
    + [nf_pipes[n]]
    + sp_dnn_l1sp_segment(n),
    'nfsp_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local f = import 'pgrapher/experiment/pdhd/funcs.jsonnet';
local fanpipe = f.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'sn_mag_nf');

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
        'gauss\\d': 'gauss',
        'wiener\\d': 'wiener',
      },
    }],
  },
}, nin=1, nout=1);

local sink = g.pnode({ type: 'DumpFrames' }, nin=1, nout=0);

// Main graph: NF + SP + DNN-ROI (+L1SP) per anode, fanned in, retagged
// (gauss%d->gauss, wiener%d->wiener) and saved as the POST-DNN result.
local graph = g.pipeline([wcls_input.adc_digits, fanpipe, retagger, wcls_output.sp_signals, sink]);

// ── Optional second output: traditional SP taken BEFORE DNN-ROI ──────
//
// An incomplete sink-like subgraph that fans in the per-anode pre-DNN
// SP frames, merges their gauss%d/wiener%d trace tags, and writes them
// through wclsFrameSaver:tradspsaver.  Its open input ports are filled
// by g.splice() below, one per broken OmnibusSigProc output edge.
//
// TagSelector is inserted before the FrameFanin so that only the
// traditional SP tags gauss%d and wiener%d are allowed into the
// tradsp save branch.  The OmnibusSigProc output carries 12 tagged
// trace sets (loose_lf, mp2/mp3_roi, tight_lf, decon_charge, the ROI
// debug sets, etc.); without this selector the full 28800-trace frame
// per anode is buffered through the fanin, roughly doubling peak SP
// memory (~8 GB).  Selecting just gauss%d/wiener%d at the tap drops the
// tradsp branch to ~5 GB.
local tradsp_preselect = [
  g.pnode({
    type: 'TagSelector',
    name: 'tradsp_preselect%d' % n,
    data: {
      tags: [
        'gauss%d' % n,
        'wiener%d' % n,
      ],
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local ofanin = g.pnode({
  type: 'FrameFanin',
  name: 'tradsp_outfanin',
  data: {
    multiplicity: std.length(tools.anodes),
    tag_rules: [
      {
        frame: { '.*': 'tradsp_outfanin' },
        trace: {
          ['gauss%d' % n]: ['gauss%d' % n],
          ['wiener%d' % n]: ['wiener%d' % n],
        },
      }
      for n in std.range(0, std.length(tools.anodes) - 1)
    ],
  },
}, nin=std.length(tools.anodes), nout=1);

local outretagger = g.pnode({
  type: 'Retagger',
  name: 'tradsp_spout',
  data: {
    tag_rules: [{
      frame: { '.*': 'tradsp_spretagger' },
      merge: {
        'gauss\\d': 'tradspgauss',
        'wiener\\d': 'tradspwiener',
      },
    }],
  },
}, nin=1, nout=1);

local osink = g.pnode({ type: 'DumpFrames', name: 'tradsp_outsink', data: {} }, nin=1, nout=0);

// The out-subgraph now has one open input port per anode (each a
// TagSelector) feeding the fanin, so it is assembled with g.intern
// rather than g.pipeline.  g.splice() fills tradsp_preselect's input
// ports from the broken OmnibusSigProc edges.
local outgr = g.intern(
  innodes=tradsp_preselect,
  centernodes=[
    ofanin,
    outretagger,
    wcls_output.tradsp_signals,
    osink,
  ],
  edges=
    [
      g.edge(tradsp_preselect[n], ofanin, 0, n)
      for n in std.range(0, std.length(tools.anodes) - 1)
    ]
    + [
      g.edge(ofanin, outretagger),
      g.edge(outretagger, wcls_output.tradsp_signals),
      g.edge(wcls_output.tradsp_signals, osink),
    ],
  name='tradsp_outgr'
);

// Break every OmnibusSigProc output edge (the raw pre-DNN SP frame) and
// insert a 2-way fanout: one branch continues into the existing DNN-ROI
// / L1SP chain, the other feeds the tradsp save subgraph above.
local edge_selector(e) = std.startsWith(e.tail.node, 'OmnibusSigProc:');
local fanout_factory(n, e) = { type: 'FrameFanout', name: 'tradsp_splice%d' % n, data: { multiplicity: 2 } };

local spliced_graph =
  if save_tradsp
  then g.splice(graph, outgr, edge_selector, fanout_factory)
  else graph;

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(spliced_graph),
  },
};

// Finally, the configuration sequence
g.uses(spliced_graph) + [app]
