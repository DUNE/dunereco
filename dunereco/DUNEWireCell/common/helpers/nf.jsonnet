local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local u = import "utils.jsonnet";

local default_dft = { type: 'FftwDFT' };

function(anode, fr, chndb, nsamples, tick=0.5*wc.us, rms_cuts=[], dft=default_dft)
    local single = {
        type: 'pdOneChannelNoise',
        name: u.idents(anode),
        uses: [dft],
        data: {
            noisedb: wc.tn(chndb),
            anode: wc.tn(anode),
            dft: wc.tn(dft),
            resmp: [
            ],
        },
    };
    // In principle there may be multiple filters of each of each
    // type.  Define them above and collect them by type here.
    local filters = {
        channel: [ single, ],
        group: [ ],
        status: [ ],
    };
    // return value
    pg.pnode({
        local name = u.idents(anode),
        type: 'OmnibusNoiseFilter',
        name: name,
        data: {
            // Nonzero forces the number of ticks in the waveform
            nticks: 0,
            // channel bin ranges are ignored
            // only when the channelmask is merged to `bad`
            maskmap: {sticky: "bad", ledge: "bad", noisy: "bad"},
            channel_filters: [wc.tn(f) for f in filters.channel],
            grouped_filters: [wc.tn(f) for f in filters.group],
            channel_status_filters: [wc.tn(f) for f in filters.status],
            noisedb: wc.tn(chndb),
            intraces: 'orig' + name,
            outtraces: 'raw' + name,
        },
    }, nin=1, nout=1, uses=[chndb, anode] + filters.channel + filters.group + filters.status)


