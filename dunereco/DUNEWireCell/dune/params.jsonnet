local common_params = import "pgrapher/dune/params.jsonnet";

common_params {
    
    adc : super.adc{
        resolution: std.extVar("Nbit"),
    },
    
    elec: super.elec {
        gain: std.extVar("elecGain")*wc.mV/wc.fC,
    }

}