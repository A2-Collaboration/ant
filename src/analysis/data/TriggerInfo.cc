#include "TriggerInfo.h"

using namespace ant::analysis::data;
using namespace std;

ostream& DAQError::Print(ostream& stream) const
{
    stream << "DAQError(Module ID=" << ModuleID << " ModuleIndex=" << ModuleIndex << " ErrorCode=" << ErrorCode << ")";
    return stream;
}

ostream &TriggerInfo::Print(ostream& stream) const
{
    stream << "TriggerInfo"
           << "(EventID=" << EventID
           << " CB Energy Sum=" << CBEnergySum << " MeV"
           << " Multipicity=" << ClusterMultiplicity
           << ")";
    for(auto& error : Errors) {
        stream << "\t" << error << "\n";
    }
    return stream;
}
