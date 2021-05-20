#pragma once

#include <string>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "dune/VLNets/data/structs/VarDict.h"

namespace VLN {

class VarExtractorBase
{
public:
    VarExtractorBase(
        const std::string &prefix,
        const std::vector<std::string> &scalarVars,
        const std::vector<std::string> &vectorVars
    );
    virtual ~VarExtractorBase() = default;
    virtual void extract(const art::Event &evt, VarDict &vars);

protected:
    virtual void extractVars(const art::Event &evt, VarDict &vars) = 0;

    void setScalarVar(
        VarDict &vars, const std::string &name, double value
    ) const;

    void appendToVectorVar(
        VarDict &vars, const std::string &name, double value
    ) const;

    void initScalarVars(
        VarDict &vars, const std::vector<std::string> &names
    ) const;

    void initVectorVars(
        VarDict &vars, const std::vector<std::string> &names
    ) const;

protected:
    std::string prefix;

    std::vector<std::string> scalarVars;
    std::vector<std::string> vectorVars;
};

}

