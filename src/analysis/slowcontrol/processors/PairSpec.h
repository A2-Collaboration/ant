#pragma once

#include "AcquScalerProcessor.h"

#include "expconfig/detectors/Trigger.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct Beam : Processor {

    AcquScalerScalar IonChamber;
    AcquScalerScalar PairSpecGate;

    Beam():
        IonChamber(expconfig::detector::Trigger::ScalerName::IonChamber),
        PairSpecGate(expconfig::detector::Trigger::ScalerName::PairSpecGate)
    {}

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override;

    virtual void PopQueue() override;

};



}}}} // namespace ant::analysis::slowcontrol::processor
