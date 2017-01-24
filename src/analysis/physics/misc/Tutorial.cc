#include "Tutorial.h"

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


Tutorial::Tutorial(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{

}

void Tutorial::ProcessEvent(const TEvent& event, manager_t& manager)
{

}

// use the classes name to register the physics class inside Ant
// this is black macro magic what's used here...but it works :)
AUTO_REGISTER_PHYSICS(Tutorial)