#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class PID_PhiAngle : public Physics {

protected:

    TH2D* pid_cb_phi_corr;
    const ant::interval<double> theta_range;

public:

    PID_PhiAngle(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event& event) override;
    void ShowResult() override;
};

}}} // namespace ant::analysis::physics