#pragma once

#include "VarExtractorBase.h"

namespace VLN {

class VLNEnergyVarExtractor : public VarExtractorBase
{
public:
    VLNEnergyVarExtractor(
        const std::string &prefix, const std::string &labelVLNEnergy
    );
    ~VLNEnergyVarExtractor() = default;

protected:
    void extractVars(const art::Event &evt, VarDict &vars) override;

private:
    std::string labelVLNEnergy;
};

}

