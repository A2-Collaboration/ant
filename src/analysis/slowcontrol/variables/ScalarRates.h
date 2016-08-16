#pragma once

#include "Variable.h"

#include "base/Detector_t.h"

#include <memory>

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace variable {

struct FreeRates : Variable {

    virtual std::list<ProcessorPtr> GetNeededProcessors() const override;

    double GetPbGlass() const;
};

}}}} // namespace ant::analysis::slowcontrol::processor
