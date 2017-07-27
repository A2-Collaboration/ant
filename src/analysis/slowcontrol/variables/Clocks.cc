#include "Clocks.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;


list<Variable::ProcessorPtr> Clocks::GetNeededProcessors() const
{
    return {Processors::Beampolmon,Processors::ExpTrigger};
}


double Clocks::GetExpClock() const
{
    if(!slowcontrol_provided)
        return 1.0;

    return Processors::ExpTrigger->Reference_1MHz.Get();
}

double Clocks::GetBeampolmonClock() const
{
    if(!slowcontrol_provided)
        return 1.0;

    return Processors::Beampolmon->Reference_1MHz.Get();
}
