#pragma once

#include "VarExtractorBase.h"

namespace VLN {

class EventAddrVarExtractor : public VarExtractorBase
{
public:
    explicit EventAddrVarExtractor(const std::string &prefix);
    ~EventAddrVarExtractor() = default;

protected:
    void extractVars(const art::Event &evt, VarDict &vars) override;
};

}

