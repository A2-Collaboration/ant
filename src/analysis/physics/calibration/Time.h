#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class Time : public Physics {

protected:

    TH2D* hTime;
    TH2D* hTimeToTagger;

    const Detector_t::Type_t DetectorType;
    bool isTagger;

public:

    Time(const Detector_t::Type_t& detectorType,
         const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event& event) override;
    void ShowResult() override;
};

}}} // namespace ant::analysis::physics