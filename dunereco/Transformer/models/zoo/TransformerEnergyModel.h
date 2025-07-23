#pragma once

#include <string>
#include <vector>

#include "dunereco/Transformer/models/pt_model/PTModel.h"
#include "dunereco/Transformer/data/structs/TransformerEnergy.h"
#include "dunereco/Transformer/data/structs/VarDict.h"

namespace transformer
{

class TransformerEnergyModel
{
private:
    PTModel model;

public:
    explicit TransformerEnergyModel(const std::string &savedir);

    TransformerEnergy predict(const VarDict &varDict);
};

}
