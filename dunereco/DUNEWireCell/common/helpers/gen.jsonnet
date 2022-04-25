// Helpers to make configuration objects for components in the "gen"
// package.  See also aux.jsonnet.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local u = import "utils.jsonnet";
local aux = import "aux.jsonnet";


{
    default_seeds: [0, 1, 2, 3, 4],

    /// Configure a random generator
    random(seeds = $.default_seeds, generator="default") :: {
        type: "Random",
        name: generator + "-" +
              std.join('-', std.map(std.toString,seeds)),
        data: {
            generator: generator,
            seeds: seeds,
        }
    },

    /// Source of depos from ideal lines of ionization
    track_depos(tracks, name="", step=1.0*wc.mm) :: pg.pnode({
        type: "TrackDepos",
        name: name,
        data: {
            step_size: step,
            tracks: tracks,
        }
    }, nin=0, nout=1),


    /// Make a drifter for all volumes.  If loadby="set" we use the
    /// depo set drifter and call this function again recursively to
    /// make the "singular" drifter.  Note: singular depos can lead to
    /// a pathological slowdown in Pgrapher execution (TbbFlow is
    /// okay).
    drifter(vols, lar, rnd=$.random(), time_offset=0,
            fluctuate=true, loadby="set") ::
        if loadby=="set"
        then pg.pnode({
            type: 'DepoSetDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1, uses=[$.drifter(vols, lar, rnd,
                                          time_offset, fluctuate,
                                          "singular")])
        else pg.pnode({
            local xregions = wc.unique_list(std.flattenArrays([
                v.faces for v in vols])),
            
            type: "Drifter",
            data: lar {
                rng: wc.tn(rnd),
                xregions: xregions,
                time_offset: time_offset,
                fluctuate: fluctuate,
            },
        }, nin=1, nout=1, uses=[rnd]),
    
    // Collect individual depos into a set.  See commen about singular
    // depos.
    bagger(daq) :: pg.pnode({
        type:'DepoBagger',
        data: {
            gate: [0, daq.nticks*daq.tick],
        },
    }, nin=1, nout=1),
    
    
    // see aux.jsonnet for individual responses

    // Return "PlaneImpactResponse" objects for all planes.
    // fr is a field response object (see fr() above).
    // srs is list of "short response" config objects, eg cer()
    // lrs is list of "long response" config objects, eg rc()
    pirs(fr, srs, lrs, dft=aux.dft) :: [ {
        type: "PlaneImpactResponse",
        name : std.toString(plane),
        data : {
            plane: plane,
            field_response: wc.tn(fr),
            short_responses: [wc.tn(r) for r in srs],
            // this needs to be big enough for convolving FR*CE
            overall_short_padding: 200*wc.us,
            long_responses: [wc.tn(r) for r in lrs],
            // this needs to be big enough to convolve RC
            long_padding: 1.5*wc.ms,
            dft: wc.tn(dft),
        },
        uses: [dft, fr] + srs + lrs,
    } for plane in [0,1,2]],

    // signal simulation
    signal(anode, pirs, daq, lar, rnd=$.random(), dft=aux.dft) ::
        pg.pipeline([
            pg.pnode({
                type:'DepoTransform',
                name: u.idents(anode),
                data: {
                    rng: wc.tn(rnd),
                    dft: wc.tn(dft),
                    anode: wc.tn(anode),
                    pirs: [wc.tn(p) for p in pirs],
                    fluctuate: true,
                    drift_speed: lar.drift_speed,
                    first_frame_number: 0,
                    readout_time: daq.nticks*daq.tick, 
                    start_time: 0,
                    tick: daq.tick,
                    nsigma: 3,
                },
            }, nin=1, nout=1, uses=pirs + [anode, rnd, dft]),

            pg.pnode({
                type: 'Reframer',
                name: u.idents(anode),
                data: {
                    anode: wc.tn(anode),
                    tags: [],
                    fill: 0.0,
                    tbin: 0,
                    toffset: 0,
                    nticks: daq.nticks,
                },
            }, nin=1, nout=1, uses=[anode])]),


    // Return a frame filter config that will add in noise.
    noise(anode, filename, daq, chstat=null, rnd=$.random(), dft=aux.dft) ::
        local cs = if std.type(chstat) == "null"
                   then {tn:"", uses:[]}
                   else {tn:wc.tn(chstat), uses:[chstat]};
        local noise_model = {
            type: "EmpiricalNoiseModel",
            name: u.idents(anode),
            data: {
                anode: wc.tn(anode),
                chanstat: cs.tn,
                spectra_file: filename,
                nsamples: daq.nticks,
                period: daq.tick,
                wire_length_scale: 1.0*wc.cm, // optimization binning
                dft: wc.tn(dft),
            },
            uses: [anode, dft] + cs.uses,
        };

        pg.pnode({
            type: "AddNoise",
            name: noise_model.name,
            data: {
                rng: wc.tn(rnd),
                model: wc.tn(noise_model),
                nsamples: daq.nticks,
                replacement_percentage: 0.02, // random optimization
                dft: wc.tn(dft),
            }}, nin=1, nout=1, uses=[rnd, noise_model, dft]),


    // digitizer simulation
    digi(anode, adc) :: 
        local apaid = anode.data.ident;
        pg.pnode({
            type: "Digitizer",
            name: '%d' % apaid,
            data : adc {
                anode: wc.tn(anode),
                frame_tag: "orig%d"%apaid,
            }
        }, nin=1, nout=1, uses=[anode]),

}
