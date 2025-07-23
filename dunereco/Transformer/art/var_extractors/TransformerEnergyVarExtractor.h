#pragma once

#include "VarExtractorBase.h"

namespace transformer {

class TransformerEnergyVarExtractor : public VarExtractorBase
{
public:
    TransformerEnergyVarExtractor(
        const std::string &prefix, const std::string &labelTransformerEnergy
    );
    ~TransformerEnergyVarExtractor() = default;

protected:
    void extractVars(const art::Event &evt, VarDict &vars) override;

private:
    std::string labelTransformerEnergy;
};

}

