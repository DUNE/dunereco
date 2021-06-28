#include "VarExtractorBase.h"
#include "utils.h"

namespace VLN {

VarExtractorBase::VarExtractorBase(
    const std::string &prefix,
    const std::vector<std::string> &scalarVars,
    const std::vector<std::string> &vectorVars
) : prefix(prefix), scalarVars(scalarVars), vectorVars(vectorVars)
{ }

void VarExtractorBase::setScalarVar(
    VarDict &vars, const std::string &name, double value
) const
{
    vars.scalar[prefix + name] = value;
}

void VarExtractorBase::appendToVectorVar(
    VarDict &vars, const std::string &name, double value
) const
{
    vars.vector[prefix + name].push_back(value);
}

void VarExtractorBase::initScalarVars(
    VarDict &vars, const std::vector<std::string> &names
) const
{
    for (auto &name : names) {
        vars.scalar[prefix + name] = -1;
    }
}

void VarExtractorBase::initVectorVars(
    VarDict &vars, const std::vector<std::string> &names
) const
{
    for (auto &name : names) {
        const std::string fullname = prefix + name;

        auto it = vars.vector.find(fullname);

        if (it != vars.vector.end()) {
            it->second.clear();
        }
        else {
            vars.vector[fullname] = {};
        }
    }
}

void VarExtractorBase::extract(const art::Event &evt, VarDict &vars)
{
    initScalarVars(vars, scalarVars);
    initVectorVars(vars, vectorVars);

    extractVars(evt, vars);
}

}

