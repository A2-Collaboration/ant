#include "ScalarRates.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;


list<Variable::ProcessorPtr> FreeRates::GetNeededProcessors() const
{
    return {Processors::Beampolmon};
}

double FreeRates::GetPbGlass() const
{
    return Processors::Beampolmon->PbGlass.Get();
}

