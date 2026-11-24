// DUNE-VD specific parameters.  This file inerets from the
// generic set of parameters and overrides things specific to PDSP.

local wc = import "wirecell.jsonnet";
local base = import "pgrapher/common/params.jsonnet";

// Bottom-drift electronics-noise file, selected by the front-end gain.  Same
// as protodunevd/params.jsonnet: the CRP strip layout is the same in PD-VD
// and in the FD-VD so the PD-VD data-driven spectra are used.  The cold
// electronics has four gain settings -- 4.7 / 7.8 / 14 / 25 mV/fC -- and
// noise-spectra files currently exist only for 7.8 and 14.  Any other gain
// aborts the configuration rather than silently using a wrong-gain spectrum.
local pdvd_bottom_noise(gain) =
    local g = gain / (wc.mV / wc.fC);
    if std.abs(g - 7.8) < 0.05 then "pdvd-bottom-noise-spectra-7d8mVfC-v1.json.bz2"
    else if std.abs(g - 14.0) < 0.05 then "pdvd-bottom-noise-spectra-14mVfC-v1.json.bz2"
    else error ("FD-VD bottom noise: no spectra file for elec.gain = " + g
                + " mV/fC.  Valid cold-electronics gain settings are"
                + " 4.7/7.8/14/25 mV/fC; spectra files exist only for 7.8 and 14.");

function(params) base {

    local ncrm = if std.objectHas(params, 'ncrm') then params.ncrm else 320,

    // CRMs [0, ncrm/2) are in the bottom drift volume and are read out by
    // the cold electronics (BDE).  CRMs [ncrm/2, ncrm) are in the top drift
    // volume and are read out by the top drift electronics (TDE).
    is_top(ident):: ident >= ncrm/2,

    // This section will be overwritten in simparams.jsonnet
    det : {

        // The current DUNE-VD goemetry has only one CRP composed by 36
        // independent CRM with side 
        // CRP is on y-z while drift is on x 
        // Only one CRP is defined in this geometry 
        // CRMs are oneside anodes     
        
        response_plane: params.response_plane,

        local bot_crp_x = -651.90*wc.cm,
        local bot_resp_x = bot_crp_x+self.response_plane,
        local bot_cathode_x = bot_crp_x+650.06*wc.cm,
        local bot_face = {
            anode:    bot_crp_x,
            response: bot_resp_x,
            cathode:  bot_cathode_x},

        local top_crp_x = 652.15*wc.cm,
        local top_resp_x = top_crp_x-self.response_plane,
        local top_cathode_x = top_crp_x - 650.06*wc.cm,
        local top_face = {
            anode:    top_crp_x,
            response: top_resp_x,
            cathode:  top_cathode_x},
       
        volumes: [
            {
                wires: n,       // anode number
                name: "crm%d"%n,
                faces: if n < ncrm/2 then [ bot_face, bot_face ]
                    else [ top_face, top_face ],
            } for n in std.range(0, ncrm-1)], // std.range is inclusive, i.e. [0, crm-1],
    },

    daq: super.daq {

        tick: 0.5*wc.us, // check this in the TDR, LArSoft
        nticks: params.nticks, //9375, // 1.6 mm/us per 0.5 us assuming 6000 mm drift leght. 

        //readout_time: self.tick*self.nticks,
        //nreadouts: 1,
        //start_time: 0.0*wc.s,
        //stop_time: self.start_time + self.nreadouts*self.readout_time,
        //first_frame_number: 0,
    },

    // ADC as in protodunevd.  These are the bottom drift values.  The top
    // drift values are applied per-CRM in the digitizer (sim.jsonnet) and
    // in sigproc (sp.jsonnet) via top_adc below.
    adc: super.adc {
        resolution: if std.objectHas(params, 'adc_resolution') then params.adc_resolution else 14,
        baselines: [1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],
        fullscale: [0.2*wc.volt, 1.6*wc.volt],
    },
    // top drift: ~1 volt baselines, 0-2volt full scale
    top_adc: self.adc {
        baselines: [1.0*wc.volt, 1.0*wc.volt, 1.0*wc.volt],
        fullscale: [0.0*wc.volt, 2.0*wc.volt],
    },

    // Electronics response as in protodunevd.  One per field response, in
    // the same order as files.fields: [bottom, top].
    local elec_gain = if std.objectHas(params, 'elecGain') then params.elecGain else 7.8*wc.mV/wc.fC,
    elecs: [
      super.elec { // bottom drift, cold electronics (BDE)
        type: "ColdElecResponse",
        gain: elec_gain,
        shaping: 2.2 * wc.us,
        postgain: 1.0,
      },
      super.elec { // top drift electronics (TDE)
        type: "JsonElecResponse",
        filename: "dunevd-coldbox-elecresp-top-psnorm_400.json.bz2",
        gain: elec_gain, // not used by JsonElecResponse
        postgain: 1.36, // 11mV/fC, 1.94 -> 14mV/fC
      },
    ],
    elec: $.elecs[0], // nominal

    // As in protodunevd: no RC response in the simulation.
    rc_resp: {
        width: 1.1*wc.ms,
        rc_layers: 0,
    },

    sim: super.sim {

        // For running in LArSoft, the simulation must be in fixed time mode. 
        fixed: true,

        // The "absolute" time (ie, relative to trigger time?) that the lower edge
        // of final readout tick #0 should correspond to.  This is a
        // "fixed" notion.
        local tick0_time = 0*wc.us,

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

    overall_short_padding: 0.2*wc.ms,
    sys_status: false,
    sys_resp: {
        start: 0.0 * wc.us,
        magnitude: 1.0,
        time_smear: 1.0 * wc.us,
    },

    files: {

        // Standard wire geometry with 2 wire planes and third dummy induction
        wires: "dunevd10kt-1x6x6-3view30deg-wires-v1.json.bz2",

        // Same CRP strip geometry as PD-VD so use its field response, one
        // per drift volume: [bottom, top].  The response_plane parameter
        // (18.1 cm) MUST match this file.
        fields: [
            "protodunevd_FR_imbalance3p_260501.json.bz2",
            "protodunevd_FR_imbalance3p_260501.json.bz2",
        ],

        // Electronics-noise spectra as in protodunevd: [bottom, top].
        // Bottom is selected by the front-end gain.
        noises: [
            pdvd_bottom_noise($.elecs[0].gain),
            "pdvd-top-noise-spectra-v3.json.bz2",
        ],

        chresp: null,

    },

}
