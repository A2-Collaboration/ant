#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class Time : public Physics {

protected:

    TH2D* hTime;
    TH2D* hTimeToF;
    TH2D* hTimeToTagger;
    TH1D* hCBTriggerTiming;

    bool isTagger;
    std::shared_ptr<Detector_t> Detector;

public:

    Time(const Detector_t::Type_t& detectorType,
         const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const TEvent& event) override;
    void ShowResult() override;
};

}}} // namespace ant::analysis::physics