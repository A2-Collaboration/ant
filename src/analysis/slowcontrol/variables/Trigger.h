#pragma once

#include "Variable.h"

#include "base/Detector_t.h"

#include <memory>

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace variable {

struct Trigger : Variable {

    virtual std::list<ProcessorPtr> GetNeededProcessors() const override;

    double GetExpLivetime() const;
    double GetExpTrigger() const;
    double GetL1Trigger() const;
};

}}}} // namespace ant::analysis::slowcontrol::processor
