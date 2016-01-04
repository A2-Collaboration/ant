#include "JustPi0.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;
using namespace std;



JustPi0::JustPi0(const string& name, PhysOptPtr opts) :
    Physics(name, opts)
{

}

void JustPi0::ProcessEvent(const Event& event)
{

}

void JustPi0::ShowResult()
{

}


AUTO_REGISTER_PHYSICS(JustPi0)
