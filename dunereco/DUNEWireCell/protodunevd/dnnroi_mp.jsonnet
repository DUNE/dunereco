// Per-anode multi-plane DNN-ROI subgraph for ProtoDUNE-VD.
//
// Uses the DNNROIFindingMultiPlane node to feed the two induction planes
// (U+V, stacked along the channel axis) into a single model call, matching
// the per-CRP stacked-induction image the local PDVD training expects
// (1, 4, ~952, 1600 -- W is dropped in training; the fully-convolutional
// MobileNetV3-UNet stretches to the toolkit's runtime channel count).
//
// The PDVD DNN-ROI models are 4-channel.  Input trace tags, in the order
// the model expects (matches the training im_tags
// frame_loose_lf / frame_mp2_roi / frame_mp3_roi / frame_gauss):
//   loose_lf{A}, mp2_roi{A}, mp3_roi{A}, gauss{A}
// NOTE: the 4th channel is gauss -- not PDHD's tight_lf.  Per-channel
// z-scale normalization is baked into the .ts, so input_scale = 1.0.
//
// Every CRP emits DNN-ROI output for both induction planes; the W
// collection plane is passed through from standard SP gauss.
//
// Output: one frame carrying two trace tags so the downstream Magnify step
// produces a standard SP-style ROOT file (hu/hv/hw_{gauss,wiener,threshold}):
//   gauss{N}  -- the DNN-ROI output (U+V from the model, W passthrough);
//                this is the "final decon" waveform.
//   wiener{N} -- the SP Wiener output, passed through unchanged WITH its
//                per-channel threshold trace_summary, so Magnify can build
//                the hu/hv/hw_threshold{N} TH1F.
// DNNROIFindingMultiPlane tags its output traces with an empty summary and
// has no summary_tag knob, so the threshold can only reach Magnify by
// carrying the SP wiener{N} traces through -- which is what the third
// branch (wiener_keep) below does.
//
// Graph: FrameFanout(3) -> [DNNROIFindingMultiPlane, W PlaneSelector,
//        wiener-keep Retagger] -> FrameFanin(3) -> Retagger (merge the three
//        per-plane DNN trace tags into gauss{N}, keep wiener{N}).

local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(anode, ts, prefix='dnnroi',
         output_scale=1.0,
         nticks=6000,
         tick_per_slice=4,
         nchunks=1,
         mask_thresh=0.5,
         nchan=4,
         debugfile='')

  local apaid = anode.data.ident;
  local prename = prefix + std.toString(apaid);

  // Input trace tags fed to the model, in the order the model expects.
  // PDVD has only 4-channel models.  Order must match the training im_tags.
  assert nchan == 4 : 'pdvd dnnroi nchan must be 4, got %d' % nchan;
  local apa_intags = [
    'loose_lf%d' % apaid,
    'mp2_roi%d' % apaid,
    'mp3_roi%d' % apaid,
    'gauss%d' % apaid,
  ];

  // The model feeds U+V (stacked) and emits DNN-ROI output for both planes.
  local mp = pg.pnode({
    type: 'DNNROIFindingMultiPlane',
    name: prename,
    data: {
      anode: wc.tn(anode),
      planes: [0, 1],  // feed model U+V (matches training)
      output_planes: [0, 1],
      intags: apa_intags,
      decon_charge_tag: 'decon_charge%d' % apaid,
      outtags: ['dnnsp%du' % apaid, 'dnnsp%dv' % apaid],
      // 4-ch PDVD models bake per-channel z-scale normalization into the
      // .ts, so they run with input_scale 1.0.
      input_scale: 1.0,
      output_scale: output_scale,
      mask_thresh: mask_thresh,
      forward: wc.tn(ts),
      nticks: nticks,
      tick_per_slice: tick_per_slice,
      nchunks: nchunks,
      debugfile: debugfile,
    },
  }, nin=1, nout=1, uses=[ts, anode]);

  // PlaneSelector that pulls gauss%d for the W plane and re-tags it.
  local w_passthrough = pg.pnode({
    type: 'PlaneSelector',
    name: prename + '_wpass',
    data: {
      anode: wc.tn(anode),
      plane: 2,
      tags: ['gauss%d' % apaid],
      tag_rules: [{
        frame: { '.*': 'DNNROIFinding' },
        trace: { ['gauss%d' % apaid]: 'dnnsp%dw' % apaid },
      }],
    },
  }, nin=1, nout=1, uses=[anode]);

  // Third branch: keep the SP Wiener traces (and their per-channel threshold
  // trace_summary) so the merged output frame can carry wiener%d alongside
  // the DNN-ROI gauss%d.  A Retagger 'trace' rule preserves the
  // trace_summary; the other SP debug tags are dropped.
  local wiener_keep = pg.pnode({
    type: 'Retagger',
    name: prename + '_wienerkeep',
    data: {
      tag_rules: [{
        trace: { ['wiener%d' % apaid]: 'wiener%d' % apaid },
      }],
    },
  }, nin=1, nout=1);

  local dnnpipes = [mp, w_passthrough, wiener_keep];
  local mult = std.length(dnnpipes);

  local dnnfanout = pg.pnode({
    type: 'FrameFanout',
    name: prename,
    data: { multiplicity: mult },
  }, nin=1, nout=mult);

  // FrameFanin merges the three branches into one frame.  It DROPS any trace
  // tag not named by a per-port 'trace' rule, so each port must list the tags
  // it brings in: port 0 the model U+V, port 1 the W passthrough, port 2 the
  // SP wiener (summary carried along).
  local dnnfanin = pg.pnode({
    type: 'FrameFanin',
    name: prename,
    data: {
      multiplicity: mult,
      tag_rules: [
        { frame: { '.*': 'dnnsp%d' % apaid },
          trace: { ['dnnsp%du' % apaid]: 'dnnsp%du' % apaid,
                   ['dnnsp%dv' % apaid]: 'dnnsp%dv' % apaid } },
        { frame: { '.*': 'dnnsp%d' % apaid },
          trace: { ['dnnsp%dw' % apaid]: 'dnnsp%dw' % apaid } },
        { frame: { '.*': 'dnnsp%d' % apaid },
          trace: { ['wiener%d' % apaid]: 'wiener%d' % apaid } },
      ],
    },
  }, nin=mult, nout=1);

  // Consolidate the three per-plane DNN-ROI trace tags into a single gauss%d
  // (the 'dnnsp%d.' regex matches dnnsp%du / dnnsp%dv / dnnsp%dw -- it assumes
  // a single-digit anode id, which holds for PDVD's 8 anodes 0-7).  wiener%d
  // is kept untouched, its threshold summary included.
  local final_retag = pg.pnode({
    type: 'Retagger',
    name: prename + '_merge',
    data: {
      tag_rules: [{
        frame: { '.*': 'dnnsp%d' % apaid },
        trace: { ['wiener%d' % apaid]: 'wiener%d' % apaid },
        merge: { ['dnnsp%d.' % apaid]: 'gauss%d' % apaid },
      }],
    },
  }, nin=1, nout=1);

  pg.intern(
    innodes=[dnnfanout],
    outnodes=[final_retag],
    centernodes=dnnpipes + [dnnfanin],
    edges=[pg.edge(dnnfanout, dnnpipes[i], i, 0) for i in std.range(0, mult - 1)]
          + [pg.edge(dnnpipes[i], dnnfanin, 0, i) for i in std.range(0, mult - 1)]
          + [pg.edge(dnnfanin, final_retag, 0, 0)],
  )
