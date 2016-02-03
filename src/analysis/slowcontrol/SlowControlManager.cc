#include "SlowControlManager.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "SlowControlVariables.h"

#include "base/Logger.h"

#include <stdexcept>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::slowcontrol;


TID SlowControlManager::min(const TID& a, const TID& b)
{
    if(a.IsInvalid())
        return b;
    if(!b.IsInvalid())
        return a;
    return std::min(a,b);
}

SlowControlManager::SlowControlManager()
{
    unsigned nRegistered = 0;
    for(VariablePtr var : Variables::All)
        /// \todo remember the processors here...
        if(!var->GetProcessors().empty())
            nRegistered++;
    LOG(INFO) << "Have " << nRegistered << " registered slowcontrol variables";
}

void SlowControlManager::ProcessEvent(const TEvent& event, physics::manager_t& manager)
{
    if(!event.Reconstructed)
        return;

    auto& reconstructed = *event.Reconstructed;

    /// \todo re-implement!
    for(const TSlowControl& sc : reconstructed.SlowControls) {
        auto entry = slowcontrol.find(sc.GetKey());
        if(entry != slowcontrol.end()) {
            buffer_t& buffer = entry->second;
            buffer.emplace(std::make_pair(reconstructed.ID, std::move(sc)));
        }
    }
}

bool SlowControlManager::isComplete() const {

    for(const auto& entry : slowcontrol) {
        if(entry.second.empty())
            return false;
    }

    return true;
}

TID SlowControlManager::GetRunUntil() const
{
    TID minimal;

    for(const auto& entry : slowcontrol) {
        const buffer_t& buffer = entry.second;
        const auto& slread = buffer.front();
        minimal = min(minimal, slread.first);
    }

    return minimal;
}

//TID SlowControlManager::UpdateSlowcontrolData(slowcontrol::SlowControl& slc)
//{

//    if(isComplete()) {
//        auto validuntil = minimal_in_buffer.IsInvalid() ? FindMinimalTID() : minimal_in_buffer;

//        minimal_in_buffer = TID();

//        for(auto& entry : slowcontrol) {
//            buffer_t& buffer = entry.second;
//            auto slread = move(buffer.front());
//            buffer.pop();

//            FillSlowControl(slc, slread.second);

//            minimal_in_buffer = min(minimal_in_buffer, slread.first);
//        }

//        return validuntil;
//    }

//    return TID();

//}
