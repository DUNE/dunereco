local common_params = import "pgrapher/common/params.jsonnet";

common_params + {
    adc : super.adc{
        resolution: 14,
    },
}