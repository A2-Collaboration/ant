#include "TaggEff.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;


list<Variable::ProcessorPtr> TaggEff::GetNeededProcessors() const
{
    return {};
}

TaggerDetector_t::taggeff_t TaggEff::Get(unsigned channel) const
{
    (void)channel;
    return {std_ext::NaN,std_ext::NaN};
}





