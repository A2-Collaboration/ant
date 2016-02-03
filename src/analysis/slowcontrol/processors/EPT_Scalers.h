#pragma once

#include "AcquScalerProcessor.h"

#include "expconfig/detectors/EPT.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct EPT_Scalers : AcquScalerVector {
    EPT_Scalers() : AcquScalerVector(expconfig::detector::EPT::ScalerName) {}
};

}}}} // namespace ant::analysis::slowcontrol::processor
