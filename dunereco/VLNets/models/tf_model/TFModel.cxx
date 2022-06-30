#include "TFModel.h"

#include <utility>

#include <boost/numeric/conversion/cast.hpp>
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

using namespace tensorflow;

template<typename T>
inline int asInt(T x)
{
    return boost::numeric_cast<int>(x);
}

TFModel::TFModel(
    const std::string &savedir,
    const std::vector<InputConfigKeys> &scalarInputKeys,
    const std::vector<InputConfigKeys> &vectorInputKeys,
    const std::vector<std::string>     &outputKeys
) : config(savedir, scalarInputKeys, vectorInputKeys, outputKeys),
    tfSession(nullptr),
    initialized(false)
{ }

Tensor TFModel::constructDummyVectorInput(
    const std::vector<std::string> &vars, float fillValue
)
{
    Tensor result(
        DT_FLOAT, TensorShape( {1, 1, asInt(vars.size())} )
    );

    auto resultData = result.tensor<float, 3>();

    for (int varIdx = 0; varIdx < asInt(vars.size()); varIdx++) {
        resultData(0, 0, varIdx) = fillValue;
    }

    return result;
}

Tensor TFModel::constructScalarInput(
    const std::unordered_map<std::string, double> &varMap,
    const std::vector<std::string>                &vars
)
{
    Tensor result(
        DT_FLOAT, TensorShape( {1, asInt(vars.size())} )
    );
    auto resultData = result.tensor<float, 2>();

    for (int varIdx = 0; varIdx < asInt(vars.size()); varIdx++) {
        resultData(0, varIdx) = varMap.at(vars[varIdx]);
    }

    return result;
}

tensorflow::Tensor TFModel::constructVectorInput(
    const std::unordered_map<std::string, std::vector<double>> &varMap,
    const std::vector<std::string>                             &vars
)
{
    if (vars.empty()) {
        return constructDummyVectorInput(vars, 0.0);
    }

    const size_t vectorSize = (varMap.empty()) ? 0 : varMap.at(vars[0]).size();

    if (vectorSize == 0) {
        /*
         * NOTE: Fake tensor with vectorSize == 1 is needed, since otherwise
         * tensorflow fails to infer graph dimensions.
         */
        return constructDummyVectorInput(vars, 0.0);
    }

    Tensor result(
        DT_FLOAT, TensorShape({ 1, asInt(vectorSize), asInt(vars.size()) })
    );
    auto resultData = result.tensor<float, 3>();

    for (int varIdx = 0; varIdx < asInt(vars.size()); varIdx++)
    {
        const auto &values = varMap.at(vars[varIdx]);

        if (values.size() != vectorSize) {
            throw std::runtime_error("Vectors have different lengths");
        }

        for (int i = 0; i < asInt(vectorSize); i++) {
            resultData(0, i, varIdx) = values[i];
        }
    }

    return result;
}

void TFModel::initTFSession() const
{
    GraphDef graph;
    Session *session = nullptr;

    /* TODO: Find NewSession version with safer memory management */
    auto status = NewSession(SessionOptions(), &session);

    if (! status.ok())
    {
        delete session;
        throw std::runtime_error(
            "Failed to initialize TF Session: " + status.ToString()
        );
    }

    tfSession.reset(session);
    session = nullptr;

    status = ReadBinaryProto(Env::Default(), config.getModelPath(), &graph);
    if (! status.ok()) {
        throw std::runtime_error(
            "Failed to load TF Graph: " + status.ToString()
        );
    }

    status = tfSession->Create(graph);
    if (! status.ok()) {
        throw std::runtime_error(
            "Failed to create TF Session: " + status.ToString()
        );
    }
}

void TFModel::ensure_initialized() const
{
    if (initialized) {
        return;
    }

    config.load();
    initTFSession();

    initialized = true;
}

std::vector<Tensor> TFModel::predict(const VarDict &vars) const
{
    ensure_initialized();

    std::vector<std::pair<std::string, Tensor>> inputs;
    std::vector<Tensor>                         outputs;

    inputs.reserve(
        config.getScalarInputs().size() + config.getVectorInputs().size()
    );

    for (const auto &inputConfig : config.getScalarInputs()) {
        inputs.emplace_back(
            inputConfig.nodeName,
            constructScalarInput(vars.scalar, inputConfig.varNames)
        );
    }

    for (const auto &inputConfig : config.getVectorInputs()) {
        inputs.emplace_back(
            inputConfig.nodeName,
            constructVectorInput(vars.vector, inputConfig.varNames)
        );
    }

    auto status = tfSession->Run(
        inputs, config.getOutputNodes(), {}, &outputs
    );

    if (! status.ok()) {
        throw std::runtime_error(
            "Failed to run TF Session: " + status.ToString()
        );
    }

    return outputs;
}

