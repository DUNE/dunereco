#include "EventRecoVarExtractor.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaHitUtils.h"

#include "utils.h"

namespace VLN {

static const std::vector<std::string> SCALAR_VARS({
    "calE", "charge", "nHits"
});

static const std::vector<std::string> VECTOR_VARS({});

EventRecoVarExtractor::EventRecoVarExtractor(
    const std::string    &prefix,
    calo::CalorimetryAlg &algCalorimetry,
    const std::string    &labelHit,
    unsigned int         plane
) : VarExtractorBase(prefix, SCALAR_VARS, VECTOR_VARS),
    algCalorimetry(algCalorimetry),
    labelHit(labelHit),
    plane(plane)
{ }

void EventRecoVarExtractor::extractVars(const art::Event &evt, VarDict &vars)
{
    using namespace dune_ana;

    const auto hits = DUNEAnaHitUtils::GetHitsOnPlane(
        DUNEAnaEventUtils::GetHits(evt, labelHit), plane
    );

    const auto chargeCalE = calcHitsChargeCalE(
        hits, evt, algCalorimetry, plane
    );

    setScalarVar(vars, "charge", chargeCalE.first);
    setScalarVar(vars, "calE",   chargeCalE.second);
    setScalarVar(vars, "nHits",  hits.size());
}

}

