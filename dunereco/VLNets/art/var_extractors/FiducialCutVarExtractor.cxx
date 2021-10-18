#include "FiducialCutVarExtractor.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace VLN {

static const std::vector<std::string> SCALAR_VARS({ "vtxContain" });
static const std::vector<std::string> VECTOR_VARS({});

FiducialCutVarExtractor::FiducialCutVarExtractor(
    const std::string &prefix,
    const fhicl::ParameterSet &pset,
    const std::string &labelGenerator
)
  : VarExtractorBase(prefix, SCALAR_VARS, VECTOR_VARS),
    labelGenerator(labelGenerator),
    containVolMaxX(pset.get<double>("ContainVolMaxX")),
    containVolMaxY(pset.get<double>("ContainVolMaxY")),
    containVolMinZ(pset.get<double>("ContainVolMinZ")),
    containVolMaxZ(pset.get<double>("ContainVolMaxZ"))
{ }

void FiducialCutVarExtractor::extractVars(const art::Event &evt, VarDict &vars)
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

    const auto &nu = mcTruth[0]->GetNeutrino();

    auto vtxX = nu.Nu().Vx();
    auto vtxY = nu.Nu().Vy();
    auto vtxZ = nu.Nu().Vz();

    setScalarVar(vars, "vtxContain",
        (
               (std::abs(vtxX) < containVolMaxX)
            && (std::abs(vtxY) < containVolMaxY)
            && (vtxZ > containVolMinZ) && (vtxZ < containVolMaxZ)
        )
    );
}

}

