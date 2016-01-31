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

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics