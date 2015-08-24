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
};

}}} // namespace ant::analysis::physics
