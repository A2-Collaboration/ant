#include "SlowcontrolManager.h"
#include "tree/TSlowControl.h"
#include "input/slowcontrol/SlowControlCreator.h"
#include <stdexcept>

using namespace ant;
using namespace ant::analysis;


TID slowcontrol::Manager::min(const TID& a, const TID& b)
{
    if(a.IsInvalid())
        return b;
    if(!b.IsInvalid())
        return a;
    return std::min(a,b);
}

void slowcontrol::Manager::SetRequiredKeys(const std::list<ant::TSlowControl::Key> keys)
{
    for(const auto& key : keys) {
        const auto entry = slowcontrol.find(key);
        if(entry == slowcontrol.end()) {
            slowcontrol[key] = buffer_t();
        }
    }
}

void slowcontrol::Manager::ProcessSlowControls(TEvent& event)
{
    /// \todo Maybe MCTrue could also have some SlowControl stuff?
    if(!event.Reconstructed)
        return;

    auto& slowcontrols = event.Reconstructed->SlowControls;

    for(TSlowControl& sc : slowcontrols) {
        auto entry = slowcontrol.find(sc.GetKey());
        if(entry != slowcontrol.end()) {
            buffer_t& buffer = entry->second;
            buffer.emplace(std::make_pair(event.Reconstructed->ID, std::move(sc)));
        }
    }

    slowcontrols.resize(0);
}

bool slowcontrol::Manager::isComplete() const {

    for(const auto& entry : slowcontrol) {
        if(entry.second.empty())
            return false;
    }

    return true;
}

TID slowcontrol::Manager::FindMinimalTID() const
{
    TID minimal;

    for(const auto& entry : slowcontrol) {
        const buffer_t& buffer = entry.second;
        const auto& slread = buffer.front();
        minimal = min(minimal, slread.first);
    }

    return minimal;
}

TID slowcontrol::Manager::UpdateSlowcontrolData(input::SlowControl& slc)
{

    if(isComplete()) {
        auto validuntil = minimal_in_buffer.IsInvalid() ? FindMinimalTID() : minimal_in_buffer;

        minimal_in_buffer = TID();

        for(auto& entry : slowcontrol) {
            buffer_t& buffer = entry.second;
            auto slread = move(buffer.front());
            buffer.pop();

            FillSlowControl(slc, slread.second);

            minimal_in_buffer = min(minimal_in_buffer, slread.first);
        }

        return validuntil;
    }

    return TID();

}
