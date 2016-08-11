#pragma once

#include "AcquScalerProcessor.h"

#include "expconfig/detectors/Trigger.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct Beampolmon : Processor {

    AcquScalerScalar Reference_1MHz;
    AcquScalerScalar PbGlass;

    Beampolmon() :
         Reference_1MHz(expconfig::detector::Trigger::ScalerName::Beampolmon_1MHz),
         PbGlass(expconfig::detector::Trigger::ScalerName::PbGlass)
    {}

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override;

    virtual void PopQueue() override;

};



}}}} // namespace ant::analysis::slowcontrol::processor
