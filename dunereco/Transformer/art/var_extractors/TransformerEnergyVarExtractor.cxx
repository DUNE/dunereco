#include "TransformerEnergyVarExtractor.h"
#include "dunereco/Transformer/data/structs/TransformerEnergy.h"

namespace transformer {

static const std::vector<std::string> SCALAR_VARS({ "primaryE", "totalE" });
static const std::vector<std::string> VECTOR_VARS({});

TransformerEnergyVarExtractor::TransformerEnergyVarExtractor(
    const std::string &prefix, const std::string &labelTransformerEnergy
)
  : VarExtractorBase(prefix, SCALAR_VARS, VECTOR_VARS),
    labelTransformerEnergy(labelTransformerEnergy)
{ }

void TransformerEnergyVarExtractor::extractVars(const art::Event &evt, VarDict &vars)
{
    auto transformerEnergy_h = evt.getHandle<TransformerEnergy>(labelTransformerEnergy);

    if (transformerEnergy_h.failedToGet()) {
        return;
    }

    setScalarVar(vars, "primaryE", transformerEnergy_h->primaryE);
    setScalarVar(vars, "totalE",   transformerEnergy_h->totalE);
}

}

