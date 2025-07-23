#pragma once

#include "TorchHandler.h"
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "dunereco/Transformer/data/structs/VarDict.h"
#include "ModelConfig.h"

namespace at {
  class Tensor;  // 👈 forward declare the dependency
}

namespace transformer{

class PTModel
{
public:
    PTModel(const std::string &savedir);
    at::Tensor predict(const VarDict &varDict);

private:
    static at::Tensor constructDummyVectorInput(
        const std::vector<std::string> &vars, float fillValue = 0.0
    );

    static at::Tensor constructScalarInput(
        const std::unordered_map<std::string, double> &varMap,
        const std::vector<std::string>                &vars
    );

    static at::Tensor constructVectorInput(
        const std::unordered_map<std::string, std::vector<double>> &varMap,
        const std::vector<std::string>                             &vars
    );

private:
    mutable ModelConfig config;
    torch::TorchHandler jitModel;
};

}