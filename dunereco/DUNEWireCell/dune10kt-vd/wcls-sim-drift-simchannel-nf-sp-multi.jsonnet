// Modal layer: "multi-sim".  The input depo set is fanned out (DepoSetFanout)
// to one finger per anode listed in process_anode_indices.  Each finger processes
// exactly one anode and TERMINATES IN ITS OWN FrameSaver(s), so each anode yields
// SEPARATE art::Event products rather than one merged product (cf. the ICARUS
// wcls-multitpc-*-refactored pattern).
//
// Per the design agreed with the user ("both per anode"): each finger emits
//   * a sim::SimChannel product (the drifted depo flux), AND
//   * a raw::RawDigit product (post-sim ADC), and
//   * a recob::Wire product (post-SP signals).
//
// The SimChannel writer sits at the HEAD of the finger (right after the per-anode
// DepoSetFilter) so it tabulates only that anode's depos over that one anode
// -> one std::vector<sim::SimChannel> product PER anode (instance "simpleSC<id>",
// converter name "postdrift<id>") rather than one merged product for all anodes.
//
// Because a wclsFrameSaver is nin=1/nout=1 (it saves a copy and passes the
// frame through), a finger is a simple linear pipe:
//
//   DepoSetFilter(anode) -> sim -> retag -> RawDigit-saver -> sp -> Wire-saver -> dump
//
// Product instance names (= WCT frame tags = art product "label"/instance):
//   * SimChannel : "simpleSC%d" % k (writer/converter name "postdrift%d")
//   * RawDigit : frame_tag_pattern % k   (e.g. "daq%d" -> "daq7")
//   * recob::Wire : "gauss%d"/"wiener%d"(/"dnnsp%d") % k -- these are the
//                   per-anode trace tags OmnibusSigProc already emits, so each
//                   anode's signals are already uniquely named.
//
// Signature: function(params, tools, fcl_params, com) -> { graph: <pnode> }

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local mf = import 'multifan.jsonnet';

function(params, tools, fcl_params, com) {

  local ftp = fcl_params.frame_tag_pattern,
  // Trick 1: when false, the per-anode RawDigit (ADC) saver is omitted so the
  //   dense ADC waveform is never put to the art::Event -- sim feeds sigproc
  //   entirely inside the WCT graph.  Only SimChannel + recob::Wire are output.
  local save_raw = fcl_params.save_rawdigits,
  // Trick 3: hydra "DeposOrBust" short-circuit (TbbFlow only -- the FrameSync
  //   merge is async).  Fingers whose filtered depo set is empty skip the whole
  //   sim+sigproc pipeline and emit an EmptyFrame instead.
  local use_hydra = fcl_params.use_hydra,

  // Build one finger for the n-th processed anode (anode ident = k).
  local finger(n) =
    local anode = tools.anodes[n];
    local k = anode.data.ident;
    local manode = com.mega_anode_for([anode]);

    // Select only this anode's depos out of the fanned-out set.
    local ds_filter = g.pnode({
      type: 'DepoSetFilter',
      name: 'ds-filter-multi-%d' % k,
      data: { anode: wc.tn(anode) },
    }, nin=1, nout=1, uses=[anode]);

    // Per-anode SimChannel writer: a depo-flux writer over THIS one anode,
    // fed the already-filtered depo set, so it emits one std::vector<SimChannel>
    // holding only this anode's channels.  It passes the depo set through to sim.
    // It lives at the HEAD of `work` (behind DeposOrBust) so a BUSTED anode skips
    // it too -> a 0-size SimChannel product, matching the RawDigit/recob::Wire.
    local scwriter = com.depoflux_writer('simpleSC%d' % k, 'postdrift%d' % k, [anode], sparse=true);

    // Optional (trick 1) RawDigit saver chain: retag the post-sim ADC frame's
    // traces to the per-anode RawDigit tag, then save raw::RawDigit.
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
    local raw_chain = if save_raw then [raw_fanin, raw_retag, raw_saver] else [];

    // Per-anode signal processing and its recob::Wire saver (sparse ROIs, trick 2).
    local sig_tags = ['gauss%d' % k, 'wiener%d' % k]
                     + (if fcl_params.use_dnnroi then ['dnnsp%d' % k] else []);
    local sp_branch = [com.sp_pipes[n]]
                      + (if fcl_params.use_dnnroi then [com.dnnroi_node(n)] else []);
    local sig_saver = com.sig_saver('spsignals%d' % k, manode, sig_tags);

    local dump = com.dumpframes('multi-dump-%d' % k);

    // The "work": SimChannel save -> sim (+ optional RawDigit save) -> sigproc ->
    // recob::Wire save.  ALL THREE art writers (SimChannel, RawDigit, recob::Wire)
    // live INSIDE this hydra-protected subgraph so a BUSTED (empty) anode skips them
    // entirely and writes ZERO-size products.  (A saver placed AFTER the FrameSync
    // merge would run on the busted EmptyFrame too, and wclsFrameSaver always emits
    // one recob::Wire per channel -> a full nchan vector of empty ROIs instead of a
    // 0-size vector, making the hydra skip invisible in the output.)  The scwriter
    // takes the depo set in and passes it through, so it leads the pipe; depos in
    // -> frame out.
    local work = g.pipeline([scwriter, com.sn_pipes[n]] + raw_chain + sp_branch + [sig_saver],
                            'multiwork%d' % k);

    // Hydra switch: DeposOrBust routes a non-empty depo set to `work` (port 0)
    // or, when empty, emits an EmptyFrame (port 1); FrameSync merges the two.
    local dorb = g.pnode({ type: 'DeposOrBust', name: 'dorb-multi-%d' % k },
                         nin=1, nout=2);
    local frame_sync = g.pnode({ type: 'FrameSync', name: 'frame-sync-multi-%d' % k },
                               nin=2, nout=1);
    local switch = g.intern(
      innodes=[dorb], outnodes=[frame_sync], centernodes=[work],
      edges=[
        g.edge(dorb, work, 0, 0),
        g.edge(work, frame_sync, 0, 0),
        g.edge(dorb, frame_sync, 1, 1),
      ], name='multiswitch%d' % k);

    local body = if use_hydra
                 then [ds_filter, switch]
                 else [ds_filter, work];

    g.pipeline(body + [dump], 'multifinger%d' % k),

  local fingers = [finger(n) for n in com.anode_iota],

  // DepoSetFanout -> per-anode self-terminating fingers.  Each finger writes its
  // own SimChannel / RawDigit / recob::Wire products (see finger() above), so
  // nothing is written before the fan-out.  The fan is built via multifan so a
  // whole-detector (>10 anode) job factors into a TREE of <=10-way DepoSetFanouts
  // -- a flat fanout would trip the TbbFlow FanoutWrapper nout<=10 assertion.
  // Only the total anode count drives the factoring (no per-layer knobs in FCL).
  local fo = mf.fanout('multi_sim_fanout', 'DepoSetFanout', std.length(fingers), max=10),
  local fanout = g.intern(
    innodes=[fo], outnodes=fingers,
    edges=[g.edge(fo, fingers[n], n, 0) for n in com.anode_iota]),

  graph: g.pipeline([
    com.wcls_input.deposet,
    com.setdrifter,
    fanout,
  ]),
}
