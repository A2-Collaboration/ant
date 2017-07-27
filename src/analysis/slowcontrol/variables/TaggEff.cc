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
    if(!slowcontrol_provided)
        return {1.0, 0.0};

    // handle non-present tagging efficiencies as MC case
    const auto& taggeffs = Processors::TaggEff->Get();
    if(taggeffs.empty())
        return {1.0, 0.0};
    return taggeffs.at(channel);
}





