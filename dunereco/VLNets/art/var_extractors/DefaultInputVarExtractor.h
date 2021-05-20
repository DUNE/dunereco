#pragma once

#include "fhiclcpp/ParameterSet.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "VarExtractorBase.h"
#include "EventAddrVarExtractor.h"
#include "EventRecoVarExtractor.h"
#include "PFParticleVarExtractor.h"

namespace VLN {

class DefaultInputVarExtractor : public VarExtractorBase
{
public:
    DefaultInputVarExtractor(
        const std::string         &prefix,
        const fhicl::ParameterSet &pset,
        unsigned int              plane = 2
    );
    ~DefaultInputVarExtractor() = default;

protected:
    void extractVars(const art::Event &evt, VarDict &vars) override;

protected:
    calo::CalorimetryAlg   algCalorimetry;

    EventAddrVarExtractor  addrVarExtractor;
    EventRecoVarExtractor  recoVarExtractor;
    PFParticleVarExtractor particleVarExtractor;
};

}

