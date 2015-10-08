#include "TriggerInfo.h"

using namespace ant::analysis::data;
using namespace std;

ostream& DAQError::Print(ostream& stream) const
{
    stream << "DAQError(Module ID=" << module_id << " ModuleIndex=" << module_index << " ErrorCode=" << error << ")";
    return stream;
}

ostream &TriggerInfo::Print(ostream& stream) const
{
    stream << "TriggerInfo"
           << "(EventId=" << event_id
           << " CB Energy Sum=" << cb_energy_sum << " MeV"
           << " Multipicity=" << cluster_multiplicity
           << ")";
    for(auto& error : errors) {
        stream << "\t" << error << "\n";
    }
    return stream;
}
