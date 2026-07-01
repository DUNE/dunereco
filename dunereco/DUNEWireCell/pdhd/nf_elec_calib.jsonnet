// PDHD noise-filter variant for ELECTRONICS-RESPONSE CALIBRATION (validation).
//
// This is a copy of pdhd/nf.jsonnet with the per-channel electronics-response
// correction (`perChannelShaper`, registered from sigproc/ResponseShaper.cxx)
// added to the OmnibusNoiseFilter channel_filters.  The shaper rescales each
// channel's spectrum by  E_nominal(f) / E_channel(f) , i.e. it divides out the
// measured per-channel cold-electronics response and re-applies the ideal one.
// After this, every channel should look as if read out by the SAME ideal cold
// electronics -- which the downstream fitter validates.
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
local sp_filters = import 'pgrapher/experiment/pdhd/sp-filters.jsonnet';

local default_dft = { type: 'FftwDFT' };

function(params, anode, tools, chndbobj, n, name='', dft=default_dft,
         debug_dump_path='', debug_dump_groups=[],
         maskmap=null) {
    local single = {
        type: 'PDHDOneChannelNoise',
        name: name,
        uses: [dft, chndbobj, anode],
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(dft),
        },
    },
    local grouped = {
        type: 'PDHDCoherentNoiseSub',
        name: name,
        uses: [dft, chndbobj, anode] + sp_filters,
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(dft),
            rms_threshold: 0.0,
            time_filters:
              if anode.data.ident == 0
              then ['Wiener_tight_U_APA1', 'Wiener_tight_V_APA1', 'Wiener_tight_W_APA1']
              else ['Wiener_tight_U', 'Wiener_tight_V', 'Wiener_tight_W'],
            lf_tighter_filter: 'ROI_tighter_lf',
            lf_loose_filter: 'ROI_loose_lf',
            debug_dump_path: debug_dump_path,
            debug_dump_groups: debug_dump_groups,
        },
    },
    local fembfilt = {
         type: 'PDHDFEMBNoiseSub',
         name: name,
         uses: [anode],
         data: {
             anode: wc.tn(anode),
             width: 50.0,
             pad_nticks: 20,
             nsigma: 3.5,
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

            maskmap: { noisy:"bad", femb_noise:"bad"},
            channel_filters: [
                // wc.tn(single),
                wc.tn(shaper),
            ],
            grouped_filters: [
            ],
            // multigroup_chanfilters: [
            //  {
            //    channelgroups: chndbobj.data.femb_negpulse_groups,
            //    filters: [wc.tn(fembfilt)],
            //  },
            // {
            //    channelgroups: chndbobj.data.groups,
            //    filters: [wc.tn(grouped)],
            //  },
            // ],

            channel_status_filters: [
            ],
            noisedb: wc.tn(chndbobj),
            intraces: '',  // '' means use all traces (wildcard)
            outtraces: 'raw%d' % n,
        },
    }, uses=[chndbobj, anode, single, grouped, fembfilt, shaper, tools.elec_resp] + pc.uses, nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
}.pipe
