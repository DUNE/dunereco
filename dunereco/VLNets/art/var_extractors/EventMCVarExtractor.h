#pragma once

#include "nusimdata/SimulationBase/MCTruth.h"
#include "VarExtractorBase.h"

namespace VLN {

class EventMCVarExtractor : public VarExtractorBase
{
public:
    explicit EventMCVarExtractor(
        const std::string &prefix,
        const std::string &labelGenerator = "generator"
    );
    ~EventMCVarExtractor() = default;

protected:
    void extractVars(const art::Event &evt, VarDict &vars) override;

private:
    std::string labelGenerator;
};

}

