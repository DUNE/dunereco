local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local sim_maker = import "pgrapher/common/sim/nodes.jsonnet";


// return some nodes, includes base sim nodes.
function(params, tools) {
    local sim = sim_maker(params, tools),

    local nanodes = std.length(tools.anodes),

    // I rue the day that we must have an (anode) X (field) cross product!
    local ductors = sim.make_detector_ductors("nominal", tools.anodes, tools.pirs[0]),


    local zippers = [sim.make_depozipper("depozipper%d"%n, tools.anodes[n], tools.pirs[0])
                     for n in std.range(0, nanodes-1)],
    local transforms = [sim.make_depotransform("depotransform%d"%n, tools.anodes[n], tools.pirs[0])
                        for n in std.range(0, nanodes-1)],
    local depos2traces = transforms,
    //local depos2traces = zippers,

    local digitizers = [
        sim.digitizer(tools.anodes[n], name="digitizer%d"%n, tag="orig%d"%n)
        for n in std.range(0,nanodes-1)],

    local reframers = [
        g.pnode({
            type: 'Reframer',
            name: 'reframer%d'%n,
            data: {
                anode: wc.tn(tools.anodes[n]),
                tags: [],           // ?? what do?
                fill: 0.0,
                tbin: params.sim.reframer.tbin,
                toffset: 0,
                nticks: params.sim.reframer.nticks,
            },
        }, nin=1, nout=1) for n in std.range(0, nanodes-1)],


    // fixme: see https://github.com/WireCell/wire-cell-gen/issues/29
    local make_noise_model = function(anode, csdb=null) {
        type: "EmpiricalNoiseModel",
        name: "empericalnoise%s"% anode.name,
        data: {
            anode: wc.tn(anode),
            dft: wc.tn(tools.dft),
            chanstat: if std.type(csdb) == "null" then "" else wc.tn(csdb),
            spectra_file: params.files.noise,
            nsamples: params.daq.nticks,
            period: params.daq.tick,
            wire_length_scale: 1.0*wc.cm, // optimization binning
        },
        uses: [anode, tools.dft] + if std.type(csdb) == "null" then [] else [csdb],
    },

    local mega_anode = {
        type: 'MegaAnodePlane',
        name: 'meganodes',
        data: {
            anodes_tn: [wc.tn(anode) for anode in tools.anodes],
        },
    },

    local add_noise = function(model, name="") g.pnode({
        type: "AddNoise",
        name: "addnoise%s"%[name],
        data: {
            rng: wc.tn(tools.random),
            dft: wc.tn(tools.dft),
            model: wc.tn(model),
            nsamples: params.daq.nticks,
            replacement_percentage: 0.02, // random optimization
        }}, nin=1, nout=1, uses=[tools.random, tools.dft, model]),

    // TODO: make this an option?
    local use_shared_model = true,
    local noises = if use_shared_model == true
    then [add_noise(make_noise_model(mega_anode), anode.name) for anode in tools.anodes]
    else [add_noise(make_noise_model(anode), anode.name) for anode in tools.anodes],
    
    local outtags = ["orig%d"%n for n in std.range(0, nanodes-1)],

    ret : {

        signal_pipelines: [g.pipeline([depos2traces[n], reframers[n],  digitizers[n]],
                                      name="simsigpipe%d"%n) for n in std.range(0, nanodes-1)],

        splusn_pipelines:  [g.pipeline([depos2traces[n], reframers[n], noises[n], digitizers[n]],
                                       name="simsignoipipe%d"%n) for n in std.range(0, nanodes-1)],
    
        signal: f.fanpipe('DepoSetFanout', self.signal_pipelines, 'FrameFanin', "simsignalgraph", outtags),
        splusn: f.fanpipe('DepoSetFanout', self.splusn_pipelines, 'FrameFanin', "simsplusngraph", outtags),

    } + sim,                    // tack on base for user sugar.
}.ret
