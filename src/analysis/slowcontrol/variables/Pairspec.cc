#include "Pairspec.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;


list<Variable::ProcessorPtr> Pairspec::GetNeededProcessors() const
{
    return {Processors::ExpTrigger,Processors::Beam};
}

double Pairspec::GetPairSpecGate() const
{
    if(!slowcontrol_provided)
        return 1.0;

    return Processors::Beam->PairSpecGate.Get() * 1.0e6 / Processors::ExpTrigger->Reference_1MHz.Get();
}
