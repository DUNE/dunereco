#include "VLNEnergyModel.h"

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

namespace VLN
{

const std::vector<InputConfigKeys> VLNEnergyModel::scalarInputKeys({
    { "input_slice", "vars_slice" }
});

const std::vector<InputConfigKeys> VLNEnergyModel::vectorInputKeys({
    { "input_png3d", "vars_png3d" }
});

const std::vector<std::string> VLNEnergyModel::outputKeys({
    "target_primary", "target_total"
});

VLNEnergyModel::VLNEnergyModel(const std::string &savedir)
    : model(savedir, scalarInputKeys, vectorInputKeys, outputKeys)
{ }

VLNEnergy VLNEnergyModel::predict(const VarDict &vars) const
{
    std::vector<tensorflow::Tensor> outputs = model.predict(vars);

    const float primaryE = outputs[0].tensor<float,2>()(0, 0);
    const float totalE   = outputs[1].tensor<float,2>()(0, 0);

    return VLNEnergy{ primaryE, totalE };
}

}
