#include "TriggerInfo.h"

std::ostream &ant::DAQError::Print(std::ostream &stream) const
{
    stream << "DAQError(Module ID=" << module_id << " ModuleIndex=" << module_index << " ErrorCode=" << error << ")";
    return stream;
}


std::ostream &ant::TriggerInfo::Print(std::ostream &stream) const
{
    stream << "TriggerInfo(CB Energy Sum=" << cb_energy_sum << " MeV Multipicity=" << cluster_multiplicity << ")";
    for(auto& error : errors) {
        stream << "\t" << error << "\n";
    }
    return stream;
}
