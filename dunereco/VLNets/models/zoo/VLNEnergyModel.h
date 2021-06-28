#pragma once

#include <string>
#include <vector>

#include "dune/VLNets/data/structs/VarDict.h"
#include "dune/VLNets/data/structs/VLNEnergy.h"
#include "dune/VLNets/models/tf_model/TFModel.h"

namespace VLN
{

class VLNEnergyModel
{
private:
    TFModel model;

    static const std::vector<InputConfigKeys> scalarInputKeys;
    static const std::vector<InputConfigKeys> vectorInputKeys;
    static const std::vector<std::string>     outputKeys;

public:
    explicit VLNEnergyModel(const std::string &savedir);

    VLNEnergy predict(const VarDict &varDict) const;
};

}
