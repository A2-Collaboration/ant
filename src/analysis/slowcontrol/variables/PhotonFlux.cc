#include "PhotonFlux.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;

list<Variable::ProcessorPtr> PhotonFlux::GetNeededProcessors()
{
    return {};
}

double PhotonFlux::Get() const
{
    return 5;
}

