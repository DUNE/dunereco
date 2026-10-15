// Modal layer: "multi-sim".  The input depo set is fanned out (DepoSetFanout)
// to one finger per TPC listed in process_tpc_indices.  Each finger processes
// exactly one TPC and TERMINATES IN ITS OWN FrameSaver(s), so each TPC yields
// SEPARATE art::Event products rather than one merged product (cf. the ICARUS
// wcls-multitpc-*-refactored pattern).
//
// Per the design agreed with the user ("both per TPC"): each finger emits BOTH
//   * a raw::RawDigit product (post-sim ADC), and
//   * a recob::Wire product (post-SP signals).
//
// Because a wclsFrameSaver is nin=1/nout=1 (it saves a copy and passes the
// frame through), a finger is a simple linear pipe:
//
//   DepoSetFilter(anode) -> sim -> retag -> RawDigit-saver -> sp -> Wire-saver -> dump
//
// Product instance names (= WCT frame tags = art product "label"/instance):
//   * RawDigit : frame_tag_pattern % tpc_id   (e.g. "daq%d" -> "daq7")
//   * recob::Wire : "gauss%d"/"wiener%d"(/"dnnsp%d") % tpc_id -- these are the
//                   per-anode trace tags OmnibusSigProc already emits, so each
//                   TPC's signals are already uniquely named.
//
// Signature: function(params, tools, fcl_params, com) -> { graph: <pnode> }

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(params, tools, fcl_params, com) {

  local ftp = fcl_params.frame_tag_pattern,

  // Build one finger for the n-th processed anode (TPC ident = k).
  local finger(n) =
    local anode = tools.anodes[n];
    local k = anode.data.ident;
    local manode = com.mega_anode_for([anode]);

    // Select only this TPC's depos out of the fanned-out set.
    local ds_filter = g.pnode({
      type: 'DepoSetFilter',
      name: 'ds-filter-multi-%d' % k,
      data: { anode: wc.tn(anode) },
    }, nin=1, nout=1, uses=[anode]);

    // Retag the post-sim ADC frame's traces to the per-TPC RawDigit tag so the
    // RawDigit saver (frame_tags=[ftp % k]) picks them up.
    local raw_fanin = g.pnode({
      type: 'FrameFanin',
      name: 'multi-raw-fanin-%d' % k,
      data: {
        multiplicity: 1,
        tags: [ftp % k],
        tag_rules: [{ frame: { '.*': ftp % k }, trace: { '.*': ftp % k } }],
      },
    }, nin=1, nout=1);
    local raw_retag = g.pnode({
      type: 'Retagger',
      name: 'multi-raw-retag-%d' % k,
      data: { tag_rules: [{ frame: { '.*': ftp % k }, merge: { '.*': ftp % k } }] },
    }, nin=1, nout=1);

    local raw_saver = com.raw_saver('simdigits%d' % k, manode, [ftp % k]);

    // Per-TPC signal processing and its recob::Wire saver.
    local sig_tags = ['gauss%d' % k, 'wiener%d' % k]
                     + (if fcl_params.use_dnnroi then ['dnnsp%d' % k] else []);
    local sp_branch = [com.sp_pipes[n]]
                      + (if fcl_params.use_dnnroi then [com.dnnroi_node(n)] else []);
    local sig_saver = com.sig_saver('spsignals%d' % k, manode, sig_tags);

    local dump = com.dumpframes('multi-dump-%d' % k);

    g.pipeline(
      [ds_filter, com.sn_pipes[n], raw_fanin, raw_retag, raw_saver]
      + sp_branch
      + [sig_saver, dump],
      'multifinger%d' % k),

  local fingers = [finger(n) for n in com.anode_iota],

  // DepoSetFanout -> per-TPC self-terminating fingers.
  local fanout = g.fan.sink('DepoSetFanout', fingers, 'multi_sim_fanout'),

  // One SimChannel product is written up front (before the fan-out).
  local depoflux = com.depoflux_writer('simpleSC'),

  graph: g.pipeline([
    com.wcls_input.deposet,
    com.setdrifter,
    depoflux,
    fanout,
  ]),
}
