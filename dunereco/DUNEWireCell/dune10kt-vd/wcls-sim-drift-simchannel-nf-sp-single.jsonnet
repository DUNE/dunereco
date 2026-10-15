// Modal layer: "single" family (single-sim / single-sp / single-sim-sp) and the
// detector-wide "full" fallback.  One combined FrameSaver per output tier is
// produced: all processed anodes are fanned out, their signal-processing frames
// fanned back in and merged (Retagger) into one product.  This is the historical
// behavior of wcls-sim-drift-simchannel-nf-sp.jsonnet, preserved verbatim here.
//
//   single-sim      : process one TPC (process_tpc_index); save RawDigit (daqtpc<i>)
//   single-sp       : process one TPC; save recob::Wire (gauss/wiener/dnnsp tpc<i>)
//   single-sim-sp   : process one TPC; save both
//   (anything else) : "full" -- process all anodes, merge into one product
//
// Signature: function(params, tools, fcl_params, com) -> { graph: <pnode> }
// where `com` is wcls-sim-drift-simchannel-nf-sp-common.jsonnet applied.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';

function(params, tools, fcl_params, com) {

  local idx = fcl_params.process_tpc_index,
  local is_single =
    fcl_params.process_mode == 'single-sim'
    || fcl_params.process_mode == 'single-sp'
    || fcl_params.process_mode == 'single-sim-sp',
  local wants_simdigits =
    fcl_params.process_mode == 'single-sim'
    || fcl_params.process_mode == 'single-sim-sp',

  local anode_iota = com.anode_iota,
  local mega_anode = com.mega_anode,

  // Combined output savers (one shared mega-anode).
  local sim_digits = com.raw_saver(
    'simdigits', mega_anode,
    if is_single then ['daqtpc%d' % idx] else ['daq']),

  local sp_signals = com.sig_saver(
    'spsignals', mega_anode,
    if is_single
    then ['gausstpc%d' % idx, 'wienertpc%d' % idx, 'dnnsptpc%d' % idx]
    else ['gauss', 'wiener', 'dnnsp']),

  // The sim (RawDigit) sub-pipe retags its single anode's traces to the
  // per-TPC "daqtpc<i>" name before the shared simdigits saver.
  local sim_retagger_fanins = [
    g.pnode({
      type: 'FrameFanin',
      name: 'sim_retagger_fanin%d' % n,
      data: {
        multiplicity: 1,
        tags: ['daqtpc%d' % n],
        tag_rules: [{ frame: { '.*': 'daqtpc%d' % n }, trace: { '.*': 'daqtpc%d' % n } }],
      },
    }, nin=1, nout=1)
    for n in [idx]
  ],
  local sim_retaggers = [
    g.pnode({
      type: 'Retagger',
      name: 'retagger-sim-%d' % n,
      data: { tag_rules: [{ frame: { '.*': 'daqtpc%d' % n }, merge: { '.*': 'daqtpc%d' % n } }] },
    }, nin=1, nout=1)
    for n in [idx]
  ],

  local full_sim_pipes = [
    g.pipeline(
      [com.sn_pipes[n]]
      + (if wants_simdigits
         then [sim_retagger_fanins[0], sim_retaggers[0], sim_digits] else []),
      'multipass%d' % n)
    for n in anode_iota
  ],

  // NOTE: the upstream config defines nf_pipes but does NOT wire the noise
  // filter into the chain (sim -> sp directly).  That (possibly-WIP) behavior is
  // preserved here verbatim; to enable NF, prepend com.nf_pipes[n] below.
  local full_sp_pipes = [
    g.pipeline(
      [com.sp_pipes[n]]
      + (if fcl_params.use_dnnroi then [com.dnnroi_node(n)] else []),
      'full_sp_pipes%d' % n)
    for n in anode_iota
  ],

  local multipass = [
    g.pipeline([full_sim_pipes[n], full_sp_pipes[n]], 'multipass%d' % n)
    for n in anode_iota
  ],

  // The hydra "switch" pipe: a DeposOrBust short-circuit so anodes with no
  // depos do not stall the flow graph.
  local make_switch_pipe = function(sim, sp, anode) {
    local d2f = g.pipeline([sim, sp]),
    local ds_filter = g.pnode({
      type: 'DepoSetFilter',
      name: 'ds-filter-switch-%d' % anode.data.ident,
      data: { anode: wc.tn(anode) },
    }, nin=1, nout=1, uses=[anode]),
    local dorb = g.pnode({
      type: 'DeposOrBust',
      name: 'dorb-switch-%d' % anode.data.ident,
    }, nin=1, nout=2),
    local frame_sync = g.pnode({
      type: 'FrameSync',
      name: 'frame-sync-switch-%d' % anode.data.ident,
    }, nin=2, nout=1),
    local ret2 = g.intern(
      innodes=[ds_filter], outnodes=[frame_sync], centernodes=[dorb, d2f],
      edges=[
        g.edge(ds_filter, dorb, 0, 0),
        g.edge(dorb, d2f, 0, 0),
        g.edge(d2f, frame_sync, 0, 0),
        g.edge(dorb, frame_sync, 1, 1),
      ]),
    local fout_bust = g.pnode({
      type: 'FrameFanout', name: 'fout-switch-bust-%d' % anode.data.ident,
      data: { multiplicity: 2 },
    }, nin=1, nout=2),
    local fout_rawdigits = g.pnode({
      type: 'FrameFanout', name: 'fout-switch-rawdigits-%d' % anode.data.ident,
      data: { multiplicity: 2 },
    }, nin=1, nout=2),
    local frame_sync_rawdigits = g.pnode({
      type: 'FrameSync', name: 'frame-sync-switch-rawdigits-%d' % anode.data.ident,
    }, nin=2, nout=1),
    local dump_rawdigits = g.pnode(
      { type: 'DumpFrames', name: 'switch-dump-rawdigits-%d' % anode.data.ident, data: {} },
      nin=1, nout=0),
    local ret3 = g.intern(
      innodes=[ds_filter], outnodes=[frame_sync],
      centernodes=[dorb, sim, sp, fout_bust, fout_rawdigits, frame_sync_rawdigits, dump_rawdigits],
      edges=[
        g.edge(ds_filter, dorb, 0, 0),
        g.edge(dorb, sim, 0, 0),
        g.edge(dorb, fout_bust, 1, 0),
        g.edge(sim, fout_rawdigits, 0, 0),
        g.edge(fout_rawdigits, sp, 0, 0),
        g.edge(sp, frame_sync, 0, 0),
        g.edge(fout_bust, frame_sync, 0, 1),
        g.edge(fout_rawdigits, frame_sync_rawdigits, 1, 0),
        g.edge(fout_bust, frame_sync_rawdigits, 1, 1),
        g.edge(frame_sync_rawdigits, dump_rawdigits, 0, 0),
      ]),
    ret: if fcl_params.save_rawdigits then ret3 else ret2,
  }.ret,

  local switch_pipes = [
    make_switch_pipe(full_sim_pipes[n], full_sp_pipes[n], tools.anodes[n])
    for n in anode_iota
  ],

  local process_pipes = if fcl_params.use_hydra then switch_pipes else multipass,

  local outtags = [],
  local tag_rules = {
    frame: { '.*': 'framefanin' },
    trace: { ['gauss%d' % anode.data.ident]: ['gauss%d' % anode.data.ident] for anode in tools.anodes }
           + { ['wiener%d' % anode.data.ident]: ['wiener%d' % anode.data.ident] for anode in tools.anodes }
           + { ['threshold%d' % anode.data.ident]: ['threshold%d' % anode.data.ident] for anode in tools.anodes }
           + { ['dnnsp%d' % anode.data.ident]: ['dnnsp%d' % anode.data.ident] for anode in tools.anodes },
  },

  local bi_manifold =
    if is_single
    then f.multifanpipe('DepoSetFanout', process_pipes, 'FrameFanin',
                        [1, 1], [1, 1], [1, 1], [1, 1], 'sn_mag', outtags, tag_rules)
    else f.multifanpipe('DepoSetFanout', process_pipes, 'FrameFanin',
                        [1, 2, 8, 32], [2, 4, 4, 10], [1, 2, 8, 32], [2, 4, 4, 10],
                        'sn_mag', outtags, tag_rules),

  local retagger = g.pnode({
    type: 'Retagger',
    data: {
      tag_rules: [{
        frame: { '.*': 'retagger' },
        merge:
          if is_single then {
            'gauss\\d+': 'gausstpc%d' % idx,
            'wiener\\d+': 'wienertpc%d' % idx,
            'dnnsp\\d+': 'dnnsptpc%d' % idx,
          } else {
            'gauss\\d+': 'gauss',
            'wiener\\d+': 'wiener',
            'dnnsp\\d+': 'dnnsp',
          },
      }],
    },
  }, nin=1, nout=1),

  local depoflux = com.depoflux_writer(
    if is_single then 'simpleSCtpc%d' % idx else 'simpleSC'),

  local sink = g.pnode({ type: 'DumpFrames', name: 'single-frame-sink', data: {} }, nin=1, nout=0),

  graph: g.pipeline([
    com.wcls_input.deposet,
    com.setdrifter,
    depoflux,
    bi_manifold,
    retagger,
    sp_signals,
    sink,
  ]),
}
