// ProtoDUNE-SP specific parameters.  This file inerets from the
// generic set of parameters and overrides things specific to PDSP.

local wc = import "wirecell.jsonnet";
local base = import "pgrapher/common/params.jsonnet";

base {
    // This section will be overwritten in simparams.jsonnet
    det : {

        // See:  wirecell-util wire-volumes protodune-wires-larsoft-v3.json.bz2
        // to help with defining these parameters.

        // between center lines
        local apa_cpa = 341.55*wc.cm,
        local cpa_thick = 50.8*wc.mm,
        local apa_w2w = 85.725*wc.mm,
        local plane_gap = 4.76*wc.mm,
        local apa_g2g = 114.3*wc.mm, 

        local apa_plane = 0.5*apa_g2g, // pick it to be at the grid wires

        // The "response" plane is where the field response functions
        // start.  Garfield calcualtions start somewhere relative to
        // something, here's where that is made concrete.  This MUST
        // match what field response functions also used.
        response_plane: 18.1*wc.cm, // relative to collection wires
                                    // synced to params.files.fields
        local res_plane = 0.5*apa_w2w + self.response_plane,

        // The cathode plane is like the anode cut off plane.  Any
        // depo not between the two is dropped prior to drifting.
        local cpa_plane = apa_cpa - 0.5*cpa_thick,

        volumes: [
            {
                // local world = 100,
                // local split = s*10, // 1: left, 2: right
                local anode = a, // physical anode number
                // wires: world + split + anode,
                // name: "anode%d"%(world + split + anode),
                wires: anode,
                name: "anode%d" %a,

                local sign = if a>3 then 1 else -1,
                local centerline = sign * apa_cpa,
                faces:
                // top drift volume
                if sign > 0
                then [
                    {
                        anode: centerline - apa_plane,
                        response: centerline - res_plane,
                        cathode: centerline - cpa_plane, 
                    },

                    {
                        anode: centerline - apa_plane,
                        response: centerline - res_plane,
                        cathode: centerline - cpa_plane, 
                    }
                ]
                // bottom drift volume
                else [
                    {
                        anode: centerline + apa_plane,
                        response: centerline + res_plane,
                        cathode: centerline + cpa_plane, 
                    },

                    {
                        anode: centerline + apa_plane,
                        response: centerline + res_plane,
                        cathode: centerline + cpa_plane, 
                    }
                ],
            } for a in std.range(0,7)
        ],

        // This describes some rough, overall bounding box.  It's not
        // directly needed but can be useful on the Jsonnet side, for
        // example when defining some simple kinematics.  It is
        // represented by a ray going from extreme corners of a
        // rectangular solid.  Again "wirecell-util wires-info" helps
        // to choose something.
        bounds : {
            tail: wc.point(-3.15, -3.42, 0, wc.m),
            head: wc.point(3.13, 3.42, 3.04, wc.m),
        }
    },

    daq: super.daq {
        nticks: 6000,
    },

    adc: super.adc {

        resolution: 14,

        // reuse ProtoDUNE SP values for bottom drift
        baselines: [1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],
        fullscale: [0.2*wc.volt, 1.6*wc.volt],
        // will rewrite top drift values in sim/digitizer and sigproc
        // FIXME: need a more elegant way
        // top drift: ~1 volt baselines, 0-2volt full scale
    },

    // This sets a relative gain at the input to the ADC.  Note, if
    // you are looking to fix SimDepoSource, you are in the wrong
    // place.  See the "scale" parameter of wcls.input.depos() defined
    // in pgrapher/common/ui/wcls/nodes.jsonnet.
    // also, see later overwriting in simparams.jsonnet
    elecs: [
      super.elec { // bottom drifter
        postgain: 1.1365, 
        shaping: 2.2 * wc.us,
      },
      super.elec { // top
        type: "JsonElecResponse",
        filename: "dunevd-coldbox-elecresp-top-psnorm_400.json.bz2",
        postgain: 1.94, // 1.94 <-> 14mV/fC
      },
    ],
    elec: $.elecs[0], // nominal 

    sim: super.sim {

        // For running in LArSoft, the simulation must be in fixed time mode. 
        fixed: true,

        // The "absolute" time (ie, in G4 time) that the lower edge of
        // of final readout tick #0 should correspond to.  This is a
        // "fixed" notion.
        local tick0_time = -250*wc.us,

        // Open the ductor's gate a bit early.
        local response_time_offset = $.det.response_plane / $.lar.drift_speed,
        local response_nticks = wc.roundToInt(response_time_offset / $.daq.tick),

        ductor : {
            nticks: $.daq.nticks + response_nticks,
            readout_time: self.nticks * $.daq.tick,
            start_time: tick0_time - response_time_offset,
        },

        // To counter the enlarged duration of the ductor, a Reframer
        // chops off the little early, extra time.  Note, tags depend on how 
        reframer: {
            tbin: response_nticks,
            nticks: $.daq.nticks,
        }
        
    },

    files: {
        wires: "protodunevd-wires-larsoft-v3.json.bz2",

        fields: [
            // "dunevdcrp2-FR-fixcoll-adjustind.json.bz2",
            // "dunevdcrp2-FR-fixcoll-adjustind.json.bz2", // repeat for top drifter
            "protodunevd_FR_3view_speed1d55.json.bz2", // remember to sync response plane position above
            "protodunevd_FR_3view_speed1d55.json.bz2",
        ],

        noises: [
            "pdvd-bottom-noise-spectra-v1.json.bz2",
            "pdvd-top-noise-spectra-v1d3.json.bz2",
        ],

        // chresp: "protodunevd-params-channel-responses-v0.json.bz2",
        chresp: null,
    },

}

