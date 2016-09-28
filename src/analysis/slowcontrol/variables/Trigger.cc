#include "Trigger.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;


list<Variable::ProcessorPtr> Trigger::GetNeededProcessors() const
{
    return {Processors::ExpTrigger};
}


double Trigger::GetExpLivetime() const
{
    return Processors::ExpTrigger->LiveCounter.Get() * 1.0 / Processors::ExpTrigger->Reference_1MHz.Get();
}

double Trigger::GetExpTrigger() const
{
    return Processors::ExpTrigger->Trigger.Get() * 1.0e6 / Processors::ExpTrigger->Reference_1MHz.Get();
}

double Trigger::GetL1Trigger() const
{
    return Processors::ExpTrigger->L1Trigger.Get() * 1.0e6 / Processors::ExpTrigger->Reference_1MHz.Get();
}
