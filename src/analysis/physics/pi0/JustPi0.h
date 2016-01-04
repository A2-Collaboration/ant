#pragma once

#include "physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class JustPi0 : public Physics {
protected:


public:
    JustPi0(const std::string& name, PhysOptPtr opts);

    virtual void ProcessEvent(const data::Event& event);
    virtual void ShowResult();
};

}}}
