#pragma once

#include "VarExtractorBase.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

namespace transformer {

class EventRecoEVarExtractor : public VarExtractorBase
{
public:
    EventRecoEVarExtractor(
        const std::string &prefix, const std::string &labelRecoE
    );
    ~EventRecoEVarExtractor() = default;

protected:
    void extractVars(const art::Event &evt, VarDict &vars) override;

private:
    std::string labelRecoE;
};

}

