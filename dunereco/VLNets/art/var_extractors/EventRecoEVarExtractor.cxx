#include "EventRecoEVarExtractor.h"

#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

namespace VLN {

static const std::vector<std::string> SCALAR_VARS({
    "nuE", "lepE", "hadE", "longestTrackContained"
});

static const std::vector<std::string> VECTOR_VARS({});

EventRecoEVarExtractor::EventRecoEVarExtractor(
    const std::string &prefix, const std::string &labelRecoE
)
  : VarExtractorBase(prefix, SCALAR_VARS, VECTOR_VARS),
    labelRecoE(labelRecoE)
{ }

void EventRecoEVarExtractor::extractVars(const art::Event &evt, VarDict &vars)
{
    auto recoE_h = evt.getHandle<dune::EnergyRecoOutput>(labelRecoE);

    if (recoE_h.failedToGet()) {
        return;
    }

    setScalarVar(vars, "nuE",  recoE_h->fNuLorentzVector.E());
    setScalarVar(vars, "lepE", recoE_h->fLepLorentzVector.E());
    setScalarVar(vars, "hadE", recoE_h->fHadLorentzVector.E());
    setScalarVar(
        vars, "longestTrackContained", recoE_h->longestTrackContained
    );
}

}

