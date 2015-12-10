#pragma once

#include "analysis/physics/Physics.h"

namespace ant {

namespace expconfig {
namespace detector {
struct TAPS;
}}

namespace analysis {
namespace physics {

class TAPS_Energy : public Physics {

protected:
    TH2D* ggIM = nullptr;
    TH3D* ggIM_mult = nullptr;
    TH2D* timing_cuts = nullptr;
    TH2D* h_pedestals = nullptr;

    std::shared_ptr<expconfig::detector::TAPS> taps_detector;

public:

    TAPS_Energy(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event& event) override;
    void ShowResult() override;
};

}}} // namespace ant::analysis::physics