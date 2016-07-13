#include "InterpolatedPulls.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;






InterpolatedPulls::InterpolatedPulls(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{

}

void InterpolatedPulls::ProcessEvent(const TEvent& event, manager_t& manager)
{

}

void InterpolatedPulls::ShowResult()
{

}

void InterpolatedPulls::Finish()
{

}

AUTO_REGISTER_PHYSICS(InterpolatedPulls)
