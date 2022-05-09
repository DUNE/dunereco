// Helper functions to build configuration objects for "auxiliary"
// components.  For now, you may find these actually implemented in
// "gen" or "sigproc" but are used elsewhere and really should be
// moved into "aux" or another, generic plugin.

local wc = import "wirecell.jsonnet";

{
    // Default DFT uses FFTW3
    dft : { type: "FftwDFT" },

    // Configure "wire" geometry and channel map to load from file
    wires(filename) :: {
        type:"WireSchemaFile",
        name: filename,
        data: {filename: filename}
    },

    // Define array of anode configurations in terms of wires and over
    // volumes (eg, as found at params.vols)
    anodes(wires, volumes) :: [ {
        type: "AnodePlane",
        name: "%d" % vol.wires,
        data: {
            ident: vol.wires,
            wire_schema: wc.tn(wires),
            faces: vol.faces,
        },
        uses: [wires]
    } for vol in volumes],


    // Configure a field respones object by filename.
    fr(filename) :: {
        type: "FieldResponse",
        name: filename,
        data: { filename: filename }
    },

    // Configure a "short" response for cold electronics model
    cer(shaping, gain, postgain, nticks, tick=0.5*wc.us, name="") ::
        {
            type: "ColdElecResponse",
            name: "",               // some dets may have more than 1
            data: {
                shaping: shaping,
                gain: gain,
                postgain: postgain,
                nticks: nticks,
                tick: tick,
            },            
        },
    
    // Configure "long" RC response
    rc(width, nticks, tick=0.5*wc.us, name="") :: {
        type: "RCResponse",
        name: name,
        data: {
            width: width,
            nticks: nticks,
            tick: tick,
        }
    },

}
