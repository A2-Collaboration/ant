#include "TestSigmaPlusFitter.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;





TestSigmaPlusFitter::TestSigmaPlusFitter(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{

}

void TestSigmaPlusFitter::ProcessEvent(const TEvent& event, manager_t&)
{

}

void TestSigmaPlusFitter::ShowResult()
{

}

AUTO_REGISTER_PHYSICS(TestSigmaPlusFitter)