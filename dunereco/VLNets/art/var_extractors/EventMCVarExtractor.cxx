#include "EventMCVarExtractor.h"

#include "nusimdata/SimulationBase/MCTruth.h"

namespace VLN {

static const std::vector<std::string> SCALAR_VARS({
    "isCC", "pdg", "mode", "lepPdg", "nuE", "lepE", "hadE"
});

static const std::vector<std::string> VECTOR_VARS({});

EventMCVarExtractor::EventMCVarExtractor(
    const std::string &prefix, const std::string &labelGenerator
)
  : VarExtractorBase(prefix, SCALAR_VARS, VECTOR_VARS),
    labelGenerator(labelGenerator)
{ }

void EventMCVarExtractor::extractVars(const art::Event &evt, VarDict &vars)
{
    std::vector<art::Ptr<simb::MCTruth>> mcTruth;

    auto mcTruth_h = evt.getHandle<std::vector<simb::MCTruth>>(labelGenerator);
    if (!mcTruth_h) {
        return;
    }

    art::fill_ptr_vector(mcTruth, mcTruth_h);
    if (mcTruth.empty()) {
        return;
    }

    const auto &nuInt = mcTruth[0]->GetNeutrino();

    setScalarVar(vars, "isCC",   (nuInt.CCNC() == 0));
    setScalarVar(vars, "pdg",    nuInt.Nu().PdgCode());
    setScalarVar(vars, "mode",   nuInt.Mode());
    setScalarVar(vars, "lepPdg", nuInt.Lepton().PdgCode());

    const double nuE  = nuInt.Nu().E();
    const double lepE = nuInt.Lepton().Momentum().T();

    setScalarVar(vars, "nuE",  nuE);
    setScalarVar(vars, "lepE", lepE);
    setScalarVar(vars, "hadE", nuE - lepE);
}

}

