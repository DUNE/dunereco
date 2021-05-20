#include "EventAddrVarExtractor.h"

namespace VLN {

static const std::vector<std::string> SCALAR_VARS({
    "run", "subRun", "event"
});

static const std::vector<std::string> VECTOR_VARS({});

EventAddrVarExtractor::EventAddrVarExtractor(const std::string &prefix)
  : VarExtractorBase(prefix, SCALAR_VARS, VECTOR_VARS)
{ }

void EventAddrVarExtractor::extractVars(const art::Event &evt, VarDict &vars)
{
    setScalarVar(vars, "run",    evt.id().run());
    setScalarVar(vars, "subRun", evt.id().subRun());
    setScalarVar(vars, "event",  evt.id().event());
}

}

