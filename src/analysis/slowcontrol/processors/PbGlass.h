#pragma once

#include "AcquScalerProcessor.h"

#include "expconfig/detectors/Trigger.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct PbGlass : Processor {

    AcquScalerScalar PbGlassAcqu;

    PbGlass() :
         PbGlassAcqu(expconfig::detector::Trigger::ScalerName::PbGlass)
    {}

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override {
        return PbGlassAcqu.ProcessEventData(recon, manager);
    }
    virtual void PopQueue() override {
        PbGlassAcqu.PopQueue();
    }
    virtual void Reset() override {
        PbGlassAcqu.Reset();
    }
};



}}}} // namespace ant::analysis::slowcontrol::processor
