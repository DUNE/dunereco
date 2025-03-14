local common_params = import "pgrapher/common/params.jsonnet";
local wc = import "wirecell.jsonnet";

common_params {
    
    adc : super.adc{
        resolution: std.extVar("Nbit"),
    },
    
    elec: super.elec {
        gain: std.extVar("elecGain")*wc.mV/wc.fC,
    }

}
