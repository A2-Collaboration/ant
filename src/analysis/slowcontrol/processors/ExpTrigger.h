#pragma once

#include "AcquScalerProcessor.h"

#include "expconfig/detectors/Trigger.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct ExpTrigger : Processor {

    AcquScalerScalar Reference_1MHz;
    AcquScalerScalar LiveCounter;
    AcquScalerScalar Trigger;
    AcquScalerScalar L1Trigger;

    ExpTrigger() :
         Reference_1MHz(expconfig::detector::Trigger::ScalerName::Exptrigger_1MHz),
         LiveCounter(expconfig::detector::Trigger::ScalerName::TotalLivetime),
         Trigger(expconfig::detector::Trigger::ScalerName::ExpTrigger),
         L1Trigger(expconfig::detector::Trigger::ScalerName::L1Trigger)
    {}

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override;

    virtual void PopQueue() override;

};



}}}} // namespace ant::analysis::slowcontrol::processor
