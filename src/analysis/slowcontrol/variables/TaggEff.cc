#include "TaggEff.h"
#include "SlowControlProcessors.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;

list<Variable::ProcessorPtr> TaggEff::GetNeededProcessors() const
{
    return {Processors::TaggEff};
}

TaggerDetector_t::taggeff_t TaggEff::Get(unsigned channel) const
{
    return Processors::TaggEff->Get().at(channel);
}





