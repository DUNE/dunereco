#pragma once

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "VarExtractorBase.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"

#include "utils.h"

namespace transformer {

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

