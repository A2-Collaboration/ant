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

};


}}} // end of ant::analysis::physics