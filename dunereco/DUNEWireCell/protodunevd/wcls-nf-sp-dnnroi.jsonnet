// WC/LS DNN-ROI entry point for ProtoDUNE-VD.
//
// This is the LArSoft (WC/LS) counterpart of the pure-WCT driver
// pgrapher/experiment/protodunevd/wct-nf-sp-dnnroi.jsonnet.  It reads
// RawDigits from the art::Event (wclsRawFrameSource), runs NF + SP +
// DNN-ROI + L1SP-after-DNN per anode, fans the per-anode frames back
// together and saves the result via wclsFrameSaver:spsaver.
//
// Built from protodunevd/wcls-nf-sp.jsonnet (the canonical WC/LS NF-SP
// entry) with the DNN-ROI + L1SP-after-DNN subgraph from
// wct-nf-sp-dnnroi.jsonnet spliced in per anode.
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

// FP32 best KD (6-ch) PDVD model.  Resolved via WIRECELL_PATH.
local dnnroi_model = 'dnnroi/pdvd/pipe_distill_transformer_6ch.ts';
local dnnroi_device = 'cpu';                  // 'cpu' or 'gpu'.  INT8 graph is CPU-only.
local dnnroi_nchan = 6;                        // PDVD 6-ch deployment
local dnnroi_concurrency = 1;
local dnnroi_nticks = 6000;
local dnnroi_tick_per_slice = 4;               // training rebin=4
local dnnroi_output_scale = 1.0;
local dnnroi_mask_thresh = 0.2;
local dnnroi_nchunks = 1;

// L1SP-after-DNN: enabled by default in 'hybrid' mode.
local use_l1sp_dnn = true;
local l1sp_pd_mode = 'hybrid';                 // 'process'|'dump'|'dnn'|'hybrid'|''
local l1sp_pd_adj_enable = true;
local l1sp_pd_adj_max_hops = 3;
local l1sp_pd_planes = null;                   // null -> envelope default [0,1] (U+V)
local l1sp_pd_dump_path = '';
local l1sp_pd_wf_dump_path = '';
local l1sp_pd_dump_all_rois = false;
// L1SP DNN tagger (mode 'dnn'||'hybrid').  Model resolved via WIRECELL_PATH.
local l1sp_pd_dnn_model = 'l1sp/pdvd/l1sp_dnn_pdvd_v1.ts';
local l1sp_pd_dnn_device = 'cpu';
local l1sp_pd_dnn_concurrency = 1;
local l1sp_pd_dnn_threshold = 0.5;             // single-threshold fallback
local l1sp_pd_dnn_threshold_bottom = 0.5;      // per-CRP override (apa < 4)
local l1sp_pd_dnn_threshold_top = 0.5;         // per-CRP override (apa >= 4)
local l1sp_pd_dnn_window_ticks = 256;
local l1sp_pd_dnn_debug_path = '';
local l1sp_pd_adj_dnn_veto = true;
// Loose-heur pre-filter overrides (DNN-chain values).
local l1sp_pd_gmax_min = 1500.0;
local l1sp_pd_min_length = 30;
local l1sp_pd_energy_frac_thr = 0.66;

local sp_maker = import 'pgrapher/experiment/protodunevd/sp.jsonnet';
// DNN-ROI input tags require OmnibusSigProc debug + multi-plane-protection.
local sp_override = { sparse: sigoutform == 'sparse' }
                    + (if use_dnnroi
                       then { use_roi_debug_mode: true, use_multi_plane_protection: true }
                       else {});

local use_resampler = std.extVar('use_resampler');

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
      tick: if use_resampler == 'true' then 512 * wc.ns else 500 * wc.ns,
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
      frame_scale: [0.005, 0.005],
      chanmaskmaps: [],
      nticks: -1,
    },
  }, nin=1, nout=1, uses=[mega_anode]),
};

local base = import 'chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: base(params, tools.anodes[n], tools.field, tools.anodes[n].data.ident) { dft: wc.tn(tools.dft) },
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in std.range(0, std.length(tools.anodes) - 1)];

local nf_maker = import 'pgrapher/experiment/protodunevd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], tools.anodes[n].data.ident, name='nf%d' % tools.anodes[n].data.ident) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp = sp_maker(params, tools, sp_override);
// L1SP OFF (l1sp_pd_mode='') inside SP so the DNN-ROI debug input tags
// survive; L1SP runs after DNN-ROI in the envelope below.
local sp_pipes = [sp.make_sigproc(a, l1sp_pd_mode='') for a in tools.anodes];

// TorchService shared by all per-anode DNN-ROI nodes.
local ts = {
  type: 'TorchService',
  name: 'dnnroi_pdvd',
  data: {
    model: dnnroi_model,
    device: dnnroi_device,
    concurrency: dnnroi_concurrency,
  },
};

local dnnroi_maker = import 'pgrapher/experiment/protodunevd/dnnroi_pp.jsonnet';
local dnnroi_inner_pipes = [dnnroi_maker(tools.anodes[n], ts,
                                         nticks=dnnroi_nticks,
                                         tick_per_slice=dnnroi_tick_per_slice,
                                         output_scale=dnnroi_output_scale,
                                         mask_thresh=dnnroi_mask_thresh,
                                         nchunks=dnnroi_nchunks,
                                         nchan=dnnroi_nchan)
                            for n in std.range(0, std.length(tools.anodes) - 1)];

// L1SP-after-DNN envelope.
local l1sp_dnn_maker = import 'pgrapher/experiment/protodunevd/l1sp_after_dnnroi.jsonnet';

// TorchService for the L1SP DNN tagger.  Built when mode is 'dnn' or
// 'hybrid'; null otherwise so the heuristic-only path doesn't pull
// libtorch a second time.
local l1sp_torch_service =
  if (l1sp_pd_mode == 'dnn' || l1sp_pd_mode == 'hybrid') then {
    type: 'TorchService',
    name: 'l1sp_dnn_pdvd',
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

local util = import 'pgrapher/experiment/protodunevd/funcs.jsonnet';
local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % tools.anodes[n].data.ident,
    data: {
      channels: util.anode_channels(tools.anodes[n].data.ident),
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
                       l1sp_pd_planes=l1sp_pd_planes,
                       l1sp_pd_dump_mode=l1sp_pd_mode,
                       l1sp_pd_dump_path=l1sp_pd_dump_path,
                       l1sp_pd_wf_dump_path=l1sp_pd_wf_dump_path,
                       l1sp_pd_dump_all_rois=l1sp_pd_dump_all_rois,
                       l1sp_pd_torch_service=l1sp_torch_service,
                       l1sp_pd_dnn_threshold=l1sp_pd_dnn_threshold,
                       l1sp_pd_dnn_threshold_bottom=l1sp_pd_dnn_threshold_bottom,
                       l1sp_pd_dnn_threshold_top=l1sp_pd_dnn_threshold_top,
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
    + (if use_resampler == 'true' && n < 4 then [resamplers[n]] else [])
    + [nf_pipes[n]]
    + sp_dnn_l1sp_segment(n),
    'nfsp_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local anode_ident = [tools.anodes[n].data.ident for n in std.range(0, std.length(tools.anodes) - 1)];
local fanin_tag_rules = [
  {
    frame: {
      '.*': 'framefanin',
    },
    trace: {
      ['gauss%d' % ind]: 'gauss%d' % ind,
      ['wiener%d' % ind]: 'wiener%d' % ind,
      ['threshold%d' % ind]: 'threshold%d' % ind,
      ['loose_lf%d' % ind]: 'loose_lf%d' % ind,
    },
  }
  for ind in anode_ident
];
local fanpipe = util.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'nfsp', [], [], fanin_tag_rules);


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
local graph = g.pipeline([wcls_input.adc_digits, fanpipe, retagger, wcls_output.sp_signals, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

// Finally, the configuration sequence
g.uses(graph) + [app]
