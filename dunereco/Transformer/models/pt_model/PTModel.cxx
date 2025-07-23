#include <torch/script.h>
#include "PTModel.h"
#include <cstdint>

namespace transformer{

PTModel::PTModel(const std::string &savedir):
    config(savedir, "scalar", "vector", "target"),
    jitModel(savedir + "/" + MODEL_NAME)
{config.load();}

torch::Tensor PTModel::constructDummyVectorInput(
    const std::vector<std::string> &vars, float fillValue
)
{
    return torch::full({1, 1, static_cast<int64_t>(vars.size())}, fillValue);
}

torch::Tensor PTModel::constructScalarInput(
    const std::unordered_map<std::string, double> &varMap,
    const std::vector<std::string>                &vars
)
{
    torch::Tensor result = torch::full({1, static_cast<int64_t>(vars.size())}, 0.0);
    for (size_t idx = 0; idx < vars.size(); idx++) {
      result.index({0, static_cast<int64_t>(idx)}) = varMap.at(vars[idx]);
    }
    return result;
}

torch::Tensor PTModel::constructVectorInput(
    const std::unordered_map<std::string, std::vector<double>> &varMap,
    const std::vector<std::string>                             &vars
)
{
    if (vars.empty()) {  // vars should never be empty!!
        return constructDummyVectorInput(vars, 0.0);
    }

    const int64_t nProngs = (varMap.empty()) ? 0 : varMap.at(vars[0]).size();

    if (nProngs == 0) {
        /*
         * NOTE: Fake tensor with nProngs == 1 is needed, since otherwise
         * pytorch fails to infer graph dimensions.
         */
        return constructDummyVectorInput(vars, 0.0);
    }

    torch::Tensor result = torch::full({1, nProngs, static_cast<int64_t>(vars.size())}, 0.0);

    for (size_t idx = 0; idx < vars.size(); idx++) {
        const auto &values = varMap.at(vars[idx]);
        if (static_cast<int64_t>(values.size()) != nProngs) {
            throw std::runtime_error("Prongs have different lengths");
        }
        for (int64_t pngIdx = 0; pngIdx < nProngs; pngIdx++) {
	  result.index({0, pngIdx, static_cast<int64_t>(idx)}) = values[pngIdx];
        }
    }

    return result;
}

torch::Tensor PTModel::predict(const VarDict &varDict)
{
    std::vector<torch::jit::IValue> inputs;
    torch::Tensor                   outputs;

    torch::Tensor vector = constructVectorInput(varDict.vector, config.getVectorInputs());
    torch::Tensor scalar = constructScalarInput(varDict.scalar, config.getScalarInputs());

    inputs.push_back(vector);
    inputs.push_back(scalar);

    outputs = jitModel.Predict(inputs).toTensor();

    return outputs;
}

}