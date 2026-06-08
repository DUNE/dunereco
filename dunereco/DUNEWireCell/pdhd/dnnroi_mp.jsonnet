// Per-anode multi-plane DNN-ROI subgraph for PDHD.
//
// Uses the new DNNROIFindingMultiPlane node to feed U+V (stacked along
// the channel axis) into a single model call, matching the input shape
// the local PDHD training expects (1, 3, 1600, 1500).
//
// Per the deployment plan in DNN_ROI_SP, APA0 emits only the U-plane
// DNN-ROI output; V and W are passed through from standard SP gauss.
// APAs 1-3 emit both U and V DNN-ROI outputs; W is passed through.
//
// Output frame trace tags per APA N:
//   - dnnspNu : DNN-ROI on U (model output)
//   - dnnspNv : DNN-ROI on V (model output for APA1-3; gauss for APA0)
//   - dnnspNw : standard SP gauss on W (always)
// And a unified frame tag "dnnspN" via the final Retagger.

local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(anode, ts, prefix='dnnroi',
         output_scale=1.0,
         nticks=6000,
         tick_per_slice=4,
         nchunks=1,
         mask_thresh=0.5,
         nchan=3,
         debugfile='')

  local apaid = anode.data.ident;
  local prename = prefix + std.toString(apaid);

  // Input trace tags fed to the model, in the order the model expects.
  // The 3-ch model (CP43.ts) uses the first three; the 6-ch KD/QAT models
  // add tight_lf/decon_charge/gauss.  Order must match the training im_tags.
  assert nchan == 3 || nchan == 6 : 'dnnroi nchan must be 3 or 6, got %d' % nchan;
  local apa_intags = [
    'loose_lf%d' % apaid,
    'mp2_roi%d' % apaid,
    'mp3_roi%d' % apaid,
    'tight_lf%d' % apaid,
    'decon_charge%d' % apaid,
    'gauss%d' % apaid,
  ][0:nchan];

  // Decide which planes get DNN-ROI outputs emitted.
  //   APA0:    output U only (V/W → standard SP gauss)
  //   APA1-3:  output U and V
  local emit_v_from_model = (apaid != 0);
  local model_output_planes = if emit_v_from_model then [0, 1] else [0];
  local model_outtags =
    if emit_v_from_model
    then ['dnnsp%du' % apaid, 'dnnsp%dv' % apaid]
    else ['dnnsp%du' % apaid];

  local mp = pg.pnode({
    type: 'DNNROIFindingMultiPlane',
    name: prename,
    data: {
      anode: wc.tn(anode),
      planes: [0, 1],  // always feed model U+V (matches training)
      output_planes: model_output_planes,
      intags: apa_intags,
      decon_charge_tag: 'decon_charge%d' % apaid,
      outtags: model_outtags,
      // 6-ch models bake per-channel z-scale normalization into the .ts,
      // so they run with input_scale 1.0; the 3-ch model needs 1/4000.
      input_scale: if nchan == 6 then 1.0 else 1.0 / 4000,
      output_scale: output_scale,
      mask_thresh: mask_thresh,
      forward: wc.tn(ts),
      nticks: nticks,
      tick_per_slice: tick_per_slice,
      nchunks: nchunks,
      debugfile: debugfile,
    },
  }, nin=1, nout=1, uses=[ts, anode]);

  // PlaneSelector that pulls gauss%d for the given plane and re-tags it.
  local pass_through(plane, outtag, name_suffix) = pg.pnode({
    type: 'PlaneSelector',
    name: prename + name_suffix,
    data: {
      anode: wc.tn(anode),
      plane: plane,
      tags: ['gauss%d' % apaid],
      tag_rules: [{
        frame: { '.*': 'DNNROIFinding' },
        trace: { ['gauss%d' % apaid]: outtag },
      }],
    },
  }, nin=1, nout=1, uses=[anode]);

  // W is always passthrough.  V is passthrough only on APA0.
  local v_passthrough = pass_through(1, 'dnnsp%dv' % apaid, '_vpass');
  local w_passthrough = pass_through(2, 'dnnsp%dw' % apaid, '_wpass');

  local sub_nodes =
    [mp]
    + (if emit_v_from_model then [] else [v_passthrough])
    + [w_passthrough];
  local mult = std.length(sub_nodes);

  local dnnfanout = pg.pnode({
    type: 'FrameFanout',
    name: prename,
    data: { multiplicity: mult },
  }, nin=1, nout=mult);

  local plane_letters =
    if emit_v_from_model then ['u', 'v', 'w']
    else ['u', 'v', 'w'];  // tags always present in output; just sourced differently

  local dnnfanin = pg.pnode({
    type: 'FrameFanin',
    name: prename,
    data: {
      multiplicity: mult,
      tag_rules: [
        { frame: { '.*': 'dnnsp%d%s' % [apaid, plane] },
          trace: { '.*': 'dnnsp%d%s' % [apaid, plane] } }
        for plane in plane_letters
      ],
    },
  }, nin=mult, nout=1);

  local retagger = pg.pnode({
    type: 'Retagger',
    name: 'dnnroi%d' % apaid,
    data: {
      tag_rules: [{
        frame: { '.*': 'dnnsp%d' % apaid },
        merge: { '.*': 'dnnsp%d' % apaid },
      }],
    },
  }, nin=1, nout=1);

  pg.intern(
    innodes=[dnnfanout],
    outnodes=[retagger],
    centernodes=sub_nodes + [dnnfanin],
    edges=[pg.edge(dnnfanout, sub_nodes[i], i, 0) for i in std.range(0, mult - 1)]
          + [pg.edge(sub_nodes[i], dnnfanin, 0, i) for i in std.range(0, mult - 1)]
          + [pg.edge(dnnfanin, retagger, 0, 0)],
  )
