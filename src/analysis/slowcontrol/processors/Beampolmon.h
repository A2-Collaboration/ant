#pragma once

#include "AcquScalerProcessor.h"

#include "expconfig/detectors/Trigger.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct Beampolmon : Processor {

    AcquScalerScalar Reference_1MHz;

    Beampolmon() :
         Reference_1MHz(expconfig::detector::Trigger::ScalerName::Beampolmon_1MHz)
    {}

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override {
        return Reference_1MHz.ProcessEventData(recon, manager);
    }
    virtual void PopQueue() override {
        Reference_1MHz.PopQueue();
    }
    virtual void Reset() override {
        Reference_1MHz.Reset();
    }

};



}}}} // namespace ant::analysis::slowcontrol::processor