#pragma once

#include "Variable.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace variable {

struct TaggerScalers : Variable {

    virtual std::list<ProcessorPtr> GetNeededProcessors() override;
};

}}}} // namespace ant::analysis::slowcontrol::processor