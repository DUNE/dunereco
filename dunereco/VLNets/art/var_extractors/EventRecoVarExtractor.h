#pragma once

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "VarExtractorBase.h"

namespace VLN {

class EventRecoVarExtractor : public VarExtractorBase
{
public:
    EventRecoVarExtractor(
        const std::string    &prefix,
        calo::CalorimetryAlg &algCalorimetry,
        const std::string    &labelHit,
        unsigned int         plane = 2
    );

    ~EventRecoVarExtractor() = default;

protected:
    void extractVars(const art::Event &evt, VarDict &vars) override;

private:
    calo::CalorimetryAlg &algCalorimetry;
    std::string          labelHit;
    unsigned int         plane;
};

}

