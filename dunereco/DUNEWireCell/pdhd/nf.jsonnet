// This provides some noise filtering related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

local default_dft = { type: 'FftwDFT' };

function(params, anode, tools, chndbobj, n, name='', dft=default_dft) {
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
        uses: [dft, chndbobj, anode],
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(dft),
            rms_threshold: 0.0,
        },
    },

//--------------------------------------------
// To add the per-channel electronics response to the workflow:

    local pc = tools.perchanresp_nameuses,
    local shaper = {
      type: 'perChannelShaper',
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
        elecresponse: wc.tn(tools.elec_resp),
      },
    },
//--------------------------------------------


    local obnf = g.pnode({
        type: 'OmnibusNoiseFilter',
        name: name,
        data: {

            // Nonzero forces the number of ticks in the waveform
            nticks: 0,

            // channel bin ranges are ignored
            // only when the channelmask is merged to `bad`
            // maskmap: {sticky: "bad", ledge: "bad", noisy: "bad"},
            channel_filters: [
                //wc.tn(single),
                wc.tn(shaper), //New
            ],
            grouped_filters: [
                wc.tn(grouped),
            ],
            channel_status_filters: [
            ],
            noisedb: wc.tn(chndbobj),
            // intraces: 'orig%d' % n,  // frame tag get all traces
            intraces: 'orig',  // frame tag get all traces
            outtraces: 'raw%d' % n,
        },
    }, uses=[chndbobj, anode, single, shaper, grouped, tools.elec_resp] + pc.uses, nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
}.pipe
