#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "dune/VLNets/data/structs/VarDict.h"
#include "ModelConfig.h"

namespace tensorflow { class Session; class Tensor; }

class TFModel
{
public:
    TFModel(
        const std::string &savedir,
        const std::vector<InputConfigKeys> &scalarInputKeys,
        const std::vector<InputConfigKeys> &vectorInputKeys,
        const std::vector<std::string>     &outputKeys
    );

    void ensure_initialized() const;
    std::vector<tensorflow::Tensor> predict(const VarDict &vars) const;

private:
    static tensorflow::Tensor constructDummyVectorInput(
        const std::vector<std::string> &vars, float fillValue = 0.0
    );

    static tensorflow::Tensor constructScalarInput(
        const std::unordered_map<std::string, double> &varMap,
        const std::vector<std::string>                &vars
    );

    static tensorflow::Tensor constructVectorInput(
        const std::unordered_map<std::string, std::vector<double>> &varMap,
        const std::vector<std::string>                             &vars
    );

    void initTFSession() const;

private:
    mutable ModelConfig config;
    mutable std::shared_ptr<tensorflow::Session> tfSession;
    mutable bool initialized;
};

