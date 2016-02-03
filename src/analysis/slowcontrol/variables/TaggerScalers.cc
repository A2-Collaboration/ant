#include "TaggerScalers.h"

#include "SlowControlProcessors.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;

list<Variable::ProcessorPtr> TaggerScalers::GetNeededProcessors()
{
    return {Processors::EPT_Scalers};
}

