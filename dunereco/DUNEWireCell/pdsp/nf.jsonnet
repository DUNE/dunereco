// This provides some noise filtering related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local gainmap = import 'pgrapher/experiment/pdsp/chndb-rel-gain.jsonnet';

local default_dft = { type: 'FftwDFT' };

function(params, anode, chndbobj, n, name='', dft=default_dft) {
    local single = {
        type: 'pdOneChannelNoise',
        name: name,
        uses: [dft, chndbobj, anode],
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(dft),
            resmp: [
                {channels: std.range(2128, 2175), sample_from: 5996},
                {channels: std.range(1520, 1559), sample_from: 5996},
                {channels: std.range( 440,  479), sample_from: 5996},
            ],
        },
    },
    local grouped = {
        type: 'mbCoherentNoiseSub',
        name: name,
        uses: [dft, chndbobj, anode],
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(dft),
            rms_threshold: 0.0,
        },
    },
    local sticky = {
        type: 'pdStickyCodeMitig',
        name: name,
        uses: [dft, chndbobj, anode],
        data: {
            extra_stky: [
                {channels: std.range(n * 2560, (n + 1) * 2560 - 1), bits: [0,1,63]},
                {channels: [4], bits: [6]  },
                {channels: [159], bits: [6]  },
                {channels: [164], bits: [36] },
                {channels: [168], bits: [7]  },
                {channels: [323], bits: [24] },
                {channels: [451], bits: [25] },
            ],
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(dft),
            stky_sig_like_val: 15.0,
            stky_sig_like_rms: 2.0,
            stky_max_len: 10,
        },
    },
    local gaincalib = {
        type: 'pdRelGainCalib',
        name: name,
        uses: [chndbobj, anode],
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            rel_gain: gainmap.rel_gain,
        },
    },

    local obnf = g.pnode({
        type: 'OmnibusNoiseFilter',
        name: name,
        data: {

            // Nonzero forces the number of ticks in the waveform
            nticks: 0,

            // channel bin ranges are ignored
            // only when the channelmask is merged to `bad`
            maskmap: {sticky: "bad", ledge: "bad", noisy: "bad"},
            channel_filters: [
                // wc.tn(sticky),
                wc.tn(single),
                // wc.tn(gaincalib),
            ],
            grouped_filters: [
                // wc.tn(grouped),
            ],
            channel_status_filters: [
            ],
            noisedb: wc.tn(chndbobj),
            intraces: 'orig%d' % n,  // frame tag get all traces
            outtraces: 'raw%d' % n,
        },
    }, uses=[chndbobj, anode, sticky, single, grouped, gaincalib], nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
}.pipe
