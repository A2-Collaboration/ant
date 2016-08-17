#include "ScalarRates.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;


list<Variable::ProcessorPtr> FreeRates::GetNeededProcessors() const
{
    return {Processors::Beampolmon,Processors::ExpTrigger};
}

double FreeRates::GetPbGlass() const
{
    return Processors::Beampolmon->PbGlass.Get() * 1.0e6 / Processors::Beampolmon->Reference_1MHz.Get();
}

double FreeRates::GetExpClock() const
{
    return Processors::ExpTrigger->Reference_1MHz.Get();
}

double FreeRates::GetBeampolmonClock() const
{
    return Processors::Beampolmon->Reference_1MHz.Get();
}

double FreeRates::GetExpLivetime() const
{
    return Processors::ExpTrigger->LiveCounter.Get() * 1.0 / Processors::ExpTrigger->Reference_1MHz.Get();
}

double FreeRates::GetExpTrigger() const
{
    return Processors::ExpTrigger->Trigger.Get() * 1.0e6 / Processors::ExpTrigger->Reference_1MHz.Get();
}

double FreeRates::GetL1Trigger() const
{
    return Processors::ExpTrigger->L1Trigger.Get() * 1.0e6 / Processors::ExpTrigger->Reference_1MHz.Get();
}
