#pragma once

#include "AcquScalerProcessor.h"

#include "expconfig/detectors/EPT.h"
#include "expconfig/detectors/Trigger.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct EPT_Scalers : AcquScalerVector {
    EPT_Scalers() : AcquScalerVector(expconfig::detector::EPT::ScalerName) {}
};

struct EPT_Or : AcquScalerScalar {
    EPT_Or() : AcquScalerScalar(expconfig::detector::Trigger::ScalerName::EPTReferenceOR){}
};

}}}} // namespace ant::analysis::slowcontrol::processor
