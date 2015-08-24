#include "DebugPhysics.h"

#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

DebugPhysics::DebugPhysics(PhysOptPtr opts): Physics("DebugPhysics", opts) {}

DebugPhysics::~DebugPhysics() {}

void DebugPhysics::ProcessEvent(const data::Event& event)
{
    VLOG(8) << event;
}

void DebugPhysics::Finish()
{
    VLOG(8) << "Nop";
}

void DebugPhysics::ShowResult()
{
    VLOG(8) << "Nop";
}

AUTO_REGISTER_PHYSICS(DebugPhysics, "DebugPhysics")
