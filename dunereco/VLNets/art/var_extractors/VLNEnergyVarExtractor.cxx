#include "VLNEnergyVarExtractor.h"
#include "dune/VLNets/data/structs/VLNEnergy.h"

namespace VLN {

static const std::vector<std::string> SCALAR_VARS({ "primaryE", "totalE" });
static const std::vector<std::string> VECTOR_VARS({});

VLNEnergyVarExtractor::VLNEnergyVarExtractor(
    const std::string &prefix, const std::string &labelVLNEnergy
)
  : VarExtractorBase(prefix, SCALAR_VARS, VECTOR_VARS),
    labelVLNEnergy(labelVLNEnergy)
{ }

void VLNEnergyVarExtractor::extractVars(const art::Event &evt, VarDict &vars)
{
    art::Handle<VLN::VLNEnergy> vlnEnergy_h;
    evt.getByLabel(labelVLNEnergy, vlnEnergy_h);

    if (vlnEnergy_h.failedToGet()) {
        return;
    }

    setScalarVar(vars, "primaryE", vlnEnergy_h->primaryE);
    setScalarVar(vars, "totalE",   vlnEnergy_h->totalE);
}

}

