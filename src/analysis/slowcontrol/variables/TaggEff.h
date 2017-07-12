#pragma once

#include "Variable.h"

#include "base/Detector_t.h"

#include <memory>

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace variable {

struct TaggEff : Variable {

    virtual std::list<ProcessorPtr> GetNeededProcessors() const override;

    TaggerDetector_t::taggeff_t Get(unsigned channel) const;
};

}}}} // namespace ant::analysis::slowcontrol::processor
