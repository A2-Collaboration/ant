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
    LOG(INFO) << event;
}

void DebugPhysics::Finish()
{
    LOG(INFO) << "Nop";
}

void DebugPhysics::ShowResult()
{
    LOG(INFO) << "Nop";
}

void DebugPhysics::Initialize(data::Slowcontrol& slowcontrol)
{
    slowcontrol.FaradayCup.Request();
}

AUTO_REGISTER_PHYSICS(DebugPhysics, "DebugPhysics")
