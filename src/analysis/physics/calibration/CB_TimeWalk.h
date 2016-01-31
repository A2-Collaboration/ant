#pragma once

#include "analysis/physics/Physics.h"

namespace ant {

namespace expconfig {
namespace detector {
struct CB;
}}

namespace analysis {
namespace physics {

class CB_TimeWalk : public Physics {

protected:
    std::shared_ptr<expconfig::detector::CB> cb_detector;
    TH3D* h_timewalk;
    TH2D* h_timewalk_overview;

public:

    CB_TimeWalk(const std::string& name, PhysOptPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics