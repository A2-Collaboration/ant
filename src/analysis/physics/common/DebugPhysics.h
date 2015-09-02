#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {


class DebugPhysics: public Physics {
public:
    DebugPhysics(PhysOptPtr opts=nullptr);
    virtual ~DebugPhysics();

    virtual void ProcessEvent(const data::Event& event);
    virtual void Finish();
    virtual void ShowResult();
    virtual void Initialize(data::Slowcontrol &slowcontrol) override;
};

class DebugPIDAlignment: public Physics {
protected:
    TH2D* angles = nullptr;
public:
    DebugPIDAlignment(PhysOptPtr opts=nullptr);
    virtual ~DebugPIDAlignment();

    virtual void ProcessEvent(const data::Event& event) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
