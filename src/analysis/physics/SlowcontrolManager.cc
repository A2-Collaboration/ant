#include "SlowcontrolManager.h"
#include "tree/TSlowControl.h"
#include "input/detail/SlowcontrolCreator.h"
#include <stdexcept>

using namespace ant;
using namespace ant::analysis;


TID slowontrol::Manager::min(const TID& a, const TID& b)
{
    if(a.IsInvalid())
        return b;
    if(!b.IsInvalid())
        return a;
    return std::min(a,b);
}

void slowontrol::Manager::SetRequiredKeys(const std::list<ant::TSlowControl::Key> keys)
{
    for(const auto& key : keys) {
        const auto entry = slowcontrol.find(key);
        if(entry == slowcontrol.end()) {
            slowcontrol[key] = buffer_t();
        }
    }
}

void slowontrol::Manager::ProcessSlowcontrol(std::unique_ptr<const TSlowControl> data)
{
    auto entry = slowcontrol.find(data->GetKey());
    if(entry != slowcontrol.end()) {
        buffer_t& buffer = entry->second;
        buffer.emplace(move(data));
    }
}

bool slowontrol::Manager::isComplete() const {

    for(const auto& entry : slowcontrol) {
        if(entry.second.empty())
            return false;
    }

    return true;
}

TID slowontrol::Manager::FindMinimalTID() const
{
    TID minimal;

    for(const auto& entry : slowcontrol) {
        const buffer_t& buffer = entry.second;
        const auto& slread = buffer.front();
        minimal = min(minimal, slread->ID);
    }

    return minimal;
}

TID slowontrol::Manager::UpdateSlowcontrolData(data::Slowcontrol& slc)
{

    if(isComplete()) {
        auto validuntil = minimal_in_buffer.IsInvalid() ? FindMinimalTID() : minimal_in_buffer;

        minimal_in_buffer = TID();

        for(auto& entry : slowcontrol) {
            buffer_t& buffer = entry.second;
            auto slread = move(buffer.front());
            buffer.pop();

            FillSlowcontrol(slc, *slread);

            minimal_in_buffer = min(minimal_in_buffer, slread->ID);
        }

        return validuntil;
    }

    return TID();

}
