#include <torch/script.h>
#include "TransformerEnergyModel.h"
#include "dunereco/Transformer/models/pt_model/PTModel.h"

namespace transformer
{

TransformerEnergyModel::TransformerEnergyModel(const std::string &savedir)
    : model(savedir)
{ }

TransformerEnergy TransformerEnergyModel::predict(const VarDict &varDict)
{
    torch::Tensor outputs = model.predict(varDict);

    const float totalE = outputs[0][0].item<float>();
    const float primaryE = outputs[0][1].item<float>();

    return TransformerEnergy{ primaryE, totalE };
}

}
