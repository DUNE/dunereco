// PDVD noise-filter variant for ELECTRONICS-RESPONSE CALIBRATION (validation).
//
// This is a copy of protodunevd/nf.jsonnet with the per-channel
// electronics-response correction (`perChannelShaper`, registered from
// sigproc/ResponseShaper.cxx) added to the OmnibusNoiseFilter channel_filters.
// The shaper rescales each channel's spectrum by  E_nominal(f) / E_channel(f) ,
// i.e. it divides out the measured per-channel cold-electronics response and
// re-applies the ideal one.  After this, every channel should look as if read
// out by the SAME ideal cold electronics -- which the downstream fitter checks.
//
// Differences from nf.jsonnet:
//   * the function takes an extra `tools` argument (needed for tools.elec_resp
//     and tools.perchanresp_nameuses).
//   * a `perChannelShaper` node is added to channel_filters.
// To activate it, the calling wcls config must set params.files.chresp to the
// fitted per-channel response JSON so tools.perchanresp -> ParamsPerChannelResponse
// is actually instantiated.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local sp_filters = import 'pgrapher/experiment/protodunevd/sp-filters.jsonnet';

local default_dft = { type: 'FftwDFT' };

function(params, anode, tools, chndbobj, n, name='', dft=default_dft,
         debug_dump_path='', debug_dump_groups=[],
         maskmap=null, shield_dump_path='',
         intraces='orig')
  {
    local single = {
      type: 'PDVDOneChannelNoise',
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
      },
    },
    // Top (_t) vs bottom (_b) anode filter suffix.  Bottom = ident 0..3,
    // top = ident 4..7.  See sp-filters.jsonnet for the registered names.
    local sfx = if anode.data.ident < 4 then '_b' else '_t',
    local grouped = {
      type: 'PDVDCoherentNoiseSub',
      name: name,
      uses: [chndbobj, anode] + sp_filters,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
        rms_threshold: 0.0,
        time_filters: ['Wiener_tight_U' + sfx,
                       'Wiener_tight_V' + sfx,
                       'Wiener_tight_W' + sfx],
        lf_tighter_filter: 'ROI_tighter_lf' + sfx,
        lf_loose_filter:   'ROI_loose_lf'   + sfx,
        debug_dump_path: debug_dump_path,
        debug_dump_groups: debug_dump_groups,
      },
    },
    local shieldcoupling_grouped = {
         type: 'PDVDShieldCouplingSub',
         name: name,
         uses: [anode],
         data: {
             anode: wc.tn(anode),
             noisedb: wc.tn(chndbobj),
             strip_length: params.files.strip_length,
             rms_threshold: 0.0,
             dump_path: shield_dump_path,
         },
     },

    // ---- Per-channel electronics-response correction (calibration) ----
    // Same singleton (ParamsPerChannelResponse) that OmnibusSigProc uses.
    local pc = tools.perchanresp_nameuses,
    local shaper = {
      type: 'perChannelShaper',
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
        elecresponse: wc.tn(tools.elec_resp),
        dft: wc.tn(dft),
      },
    },

    local obnf = g.pnode({
      type: 'OmnibusNoiseFilter',
      name: name,
      data: {

        // Nonzero forces the number of ticks in the waveform
        nticks: 0,

        maskmap: { noisy:"bad"},
        channel_filters: [
          // wc.tn(single),
          wc.tn(shaper),
        ],
            //  multigroup_chanfilters: [

            // {
            //    channelgroups: chndbobj.data.groups,
            //    filters: [wc.tn(grouped)],
            //  },

            // ]
            // // only apply to top
            // +if anode.data.ident > 3 then[
            //   {
            //     channelgroups: chndbobj.data.top_u_groups,
            //     filters: [wc.tn(shieldcoupling_grouped)],
            //   },
            // ]else []
            // ,
        grouped_filters: [
        ],
        channel_status_filters: [
        ],
        noisedb: wc.tn(chndbobj),
        intraces: intraces,  // frame tag get all traces ('' = wildcard)
        outtraces: 'raw%d' % anode.data.ident,
      },
    }, uses=[chndbobj, anode, single, grouped, shieldcoupling_grouped, shaper, tools.elec_resp] + pc.uses, nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
  }.pipe
