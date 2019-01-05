// This provides signal processing related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// BIG FAT FIXME: we are taking from uboone.  If PDSP needs tuning do
// four things: 0) read this comment, 1) cp this file into pdsp/, 2)
// fix the import and 3) delete this comment.
local spfilt = import 'pgrapher/experiment/pdsp/sp-filters.jsonnet';

function(params, tools, override = {}) {

  local pc = tools.perchanresp_nameuses,

  // pDSP needs a per-anode sigproc
  make_sigproc(anode, n, name=null):: g.pnode({
    type: 'OmnibusSigProc',
    name:
      if std.type(name) == 'null'
      then anode.name + 'sigproc%d' % n
      else name,

    data: {
      // Many parameters omitted here.
      anode: wc.tn(anode),
      field_response: wc.tn(tools.field),
      per_chan_resp: pc.name,
      fft_flag: 0,  // 1 is faster but higher memory, 0 is slightly slower but lower memory
      postgain: 1,  // default 1.2
      ADC_mV: 4096 / (1400.0 * wc.mV),  // default 4096/2000
      r_fake_signal_low_th: 400,  // default 500
      r_fake_signal_high_th: 800,  // default 1000
      r_fake_signal_low_th_ind_factor: 1.5,  // default 1
      r_fake_signal_high_th_ind_factor: 1.5,  // default 1
      troi_col_th_factor: 5.0,  // default 5
      troi_ind_th_factor: 3.5,  // default 3
      r_th_factor: 3.5,  // default 3
      wiener_tag: 'wiener%d' % n,
      wiener_threshold_tag: 'threshold%d' % n,
      gauss_tag: 'gauss%d' % n,
    } + override,
  }, nin=1, nout=1, uses=[anode, tools.field] + pc.uses + spfilt),

}
