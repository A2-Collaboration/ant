// better than using the standard ifndef include guards...
#pragma once

// physics classes need to derive from this interface
#include "physics/Physics.h"

// the physics classes reside in this nested namespace
namespace ant {
namespace analysis {
namespace physics {

// choose a nice name for your class
class Tutorial : public Physics {

public:
    // physics need to implement this public constructor...
    Tutorial(const std::string& name, OptionsPtr opts);

    // ...and the following method (use override keyword)
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

};


}}} // end of ant::analysis::physics