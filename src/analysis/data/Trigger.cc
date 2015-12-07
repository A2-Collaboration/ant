#include "Trigger.h"

using namespace ant::analysis::data;
using namespace std;

ostream& DAQError::Print(ostream& stream) const
{
    stream << "DAQError(Module ID=" << ModuleID << " ModuleIndex=" << ModuleIndex << " ErrorCode=" << ErrorCode << ")";
    return stream;
}

ostream& Trigger_t::Print(ostream& stream) const
{
    stream << "Trigger"
           << "(EventID=" << EventID
           << " CB Energy Sum=" << CBEnergySum << " MeV"
           << " Multipicity=" << ClusterMultiplicity
           << ")";
    for(auto& error : Errors) {
        stream << "\t" << error << "\n";
    }
    return stream;
}
