#pragma once

#include "AcquScalerProcessor.h"

#include "expconfig/detectors/Tagger.h"
#include "expconfig/detectors/Trigger.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct Tagger_Scalers : AcquScalerVector {
    Tagger_Scalers() : AcquScalerVector(expconfig::detector::Tagger::ScalerName) {}
};

struct Tagger_Or : AcquScalerScalar {
    Tagger_Or() : AcquScalerScalar(expconfig::detector::Trigger::ScalerName::TaggerReferenceOR){}
};

}}}} // namespace ant::analysis::slowcontrol::processor
