#include "DefaultInputVarExtractor.h"

namespace VLN {

DefaultInputVarExtractor::DefaultInputVarExtractor(
    const std::string &prefix,
    const fhicl::ParameterSet &pset,
    unsigned int plane
) : VarExtractorBase(prefix, {}, {}),
    algCalorimetry(pset.get<fhicl::ParameterSet>("AlgCalorimetry")),
    addrVarExtractor(prefix + "addr."),
    recoVarExtractor(
        prefix + "event.",
        algCalorimetry,
        pset.get<std::string>("LabelHit"),
        plane
    ),
    particleVarExtractor(
        prefix + "particle.",
        algCalorimetry,
        pset.get<std::string>("LabelPFPModule"),
        pset.get<std::string>("LabelPFPTrack"),
        pset.get<std::string>("LabelPFPShower"),
        plane
    )
{ }

void DefaultInputVarExtractor::extractVars(
    const art::Event &evt, VarDict &vars
)
{
    addrVarExtractor    .extract(evt, vars);
    recoVarExtractor    .extract(evt, vars);
    particleVarExtractor.extract(evt, vars);
}

}

