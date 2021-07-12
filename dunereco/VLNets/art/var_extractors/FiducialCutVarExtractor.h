#pragma once

#include "fhiclcpp/ParameterSet.h"
#include "VarExtractorBase.h"

namespace VLN {

class FiducialCutVarExtractor : public VarExtractorBase
{
public:
    FiducialCutVarExtractor(
        const std::string &prefix,
        const fhicl::ParameterSet &pset,
        const std::string &labelGenerator = "generator"
    );
    ~FiducialCutVarExtractor() = default;

protected:
    void extractVars(const art::Event &evt, VarDict &vars) override;

private:
    std::string labelGenerator;

    double containVolMaxX;
    double containVolMaxY;

    double containVolMinZ;
    double containVolMaxZ;
};

}

