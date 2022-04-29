// Helpers to construct config objects for components related to
// signal processing.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local aux = import "aux.jsonnet";


// Signal processing.
//
// Note, spfilt are a list of filter objects which MUST match
// hard-wired names in the C++, sorry.  See, eg
// pgrapher/experiment/pdsp/sp-filters.jsonnet.
function(anode, fieldresp, elecresp, spfilt, adcpermv, perchan=null, dft=aux.dft, override={}) 
    local apaid = anode.data.ident;

    // if perchan file name is given we need to add this to a
    // "uses" list and to set the name in OmnibusSigProc.
    local pcresp = {
        type: "PerChannelResponse",
        data: {
            filename: perchan
        }
    };
    local pc = if std.type(perchan) == "null"
               then {tn:"", uses:[]}
               else {tn:wc.tn(pcresp), uses:[pcresp]};

    // return
    pg.pnode({
        type: 'OmnibusSigProc',
        name: '%d' % apaid,
        data: {
            /**  
                *  Optimized SP parameters (May 2019)
                *  Associated tuning in sp-filters.jsonnet
                */
            anode: wc.tn(anode),
            dft: wc.tn(dft),
            field_response: wc.tn(fieldresp),
            elecresponse: wc.tn(elecresp),
            ftoffset: 0.0, // default 0.0
            ctoffset: 1.0*wc.microsecond, // default -8.0
            per_chan_resp: pc.tn,
            fft_flag: 0,  // 1 is faster but higher memory, 0 is slightly slower but lower memory
            postgain: 1.0,  // default 1.2
            ADC_mV: adcpermv, // 4096 / (1400.0 * wc.mV), 
            troi_col_th_factor: 5.0,  // default 5
            troi_ind_th_factor: 3.0,  // default 3
            lroi_rebin: 6, // default 6
            lroi_th_factor: 3.5, // default 3.5
            lroi_th_factor1: 0.7, // default 0.7
            lroi_jump_one_bin: 1, // default 0

            r_th_factor: 3.0,  // default 3
            r_fake_signal_low_th: 375,  // default 500
            r_fake_signal_high_th: 750,  // default 1000
            r_fake_signal_low_th_ind_factor: 1.0,  // default 1
            r_fake_signal_high_th_ind_factor: 1.0,  // default 1      
            r_th_peak: 3.0, // default 3.0
            r_sep_peak: 6.0, // default 6.0
            r_low_peak_sep_threshold_pre: 1200, // default 1200

            // frame tags
            wiener_tag: 'wiener%d' % apaid,
            wiener_threshold_tag: 'threshold%d' % apaid,
            decon_charge_tag: 'decon_charge%d' % apaid,
            gauss_tag: 'gauss%d' % apaid,

            use_roi_debug_mode: false,
            tight_lf_tag: 'tight_lf%d' % apaid,
            loose_lf_tag: 'loose_lf%d' % apaid,
            cleanup_roi_tag: 'cleanup_roi%d' % apaid,
            break_roi_loop1_tag: 'break_roi_1st%d' % apaid,
            break_roi_loop2_tag: 'break_roi_2nd%d' % apaid,
            shrink_roi_tag: 'shrink_roi%d' % apaid,
            extend_roi_tag: 'extend_roi%d' % apaid,

            use_multi_plane_protection: false,
            mp3_roi_tag: 'mp3_roi%d' % apaid,
            mp2_roi_tag: 'mp2_roi%d' % apaid,
            
            isWrapped: false,
            // process_planes: [0, 2],
        } + override
    }, nin=1, nout=1, uses=[anode, dft, fieldresp, elecresp] + pc.uses + spfilt)
