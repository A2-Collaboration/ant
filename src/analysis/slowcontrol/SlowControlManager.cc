#include "SlowControlManager.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "SlowControlVariables.h"

#include "base/Logger.h"

#include <stdexcept>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::slowcontrol;


void SlowControlManager::AddProcessor(std::shared_ptr<Processor> p)
{
    slowcontrol.insert({p, buffer_t()});
}

TID SlowControlManager::min(const TID& a, const TID& b)
{
    if(a.IsInvalid())
        return b;
    if(!b.IsInvalid())
        return a;
    return std::min(a,b);
}

TID SlowControlManager::PopBuffers(const TID& id)
{
    TID next;

    all_complete = true;
    for(auto& s : slowcontrol) {
        if(!s.second.empty() && s.second.front() == id) {
            s.second.pop();
            s.first->PopQueue();
        }

        if(!s.second.empty()) {
            next = min(next, s.second.front());
        } else {
            all_complete = false;
        }
    }

    return next;
}

SlowControlManager::SlowControlManager()
{
    unsigned nRegistered = 0;
    for(VariablePtr var : Variables::All)


        if(!var->GetProcessors().empty()) {

            nRegistered++;

            for(auto& p : var->GetProcessors()) {
                AddProcessor(p);
            }

        }

    LOG(INFO) << "Have " << nRegistered << " registered slowcontrol variables and " << slowcontrol.size() << " processors";
}

bool SlowControlManager::ProcessEvent(TEventPtr event)
{
    if(!event->HasReconstructed())
        return true;

    TEventData& reconstructed = event->Reconstructed();

    physics::manager_t manager;

    all_complete = true;
    bool wants_skip = false;

    for(auto& sl : slowcontrol) {

        const auto result = sl.first->ProcessEventData(reconstructed, manager);

        if(result == slowcontrol::Processor::return_t::Complete) {
            sl.second.push(reconstructed.ID);
        }
        else if(result == slowcontrol::Processor::return_t::Skip) {
            wants_skip = true;
        }

        all_complete &= !sl.second.empty();
    }

    // a skipped event could still be saved in order to trigger
    // slow control processsors (see for example AcquScalerProcessor),
    // but should NOT be processed by physics classes. Mark the event accordingly.
    if(wants_skip && !event->SavedForSlowControls)
        event->SavedForSlowControls = manager.saveEvent;

    eventbuffer.emplace(manager.saveEvent, move(event));

    return all_complete;
}

event_t SlowControlManager::PopEvent() {

    if(!all_complete || eventbuffer.empty())
        return {};

    auto i = std::move(eventbuffer.front());
    eventbuffer.pop();

    if(i.Event->Reconstructed().ID == changepoint) {
        changepoint = PopBuffers(changepoint);
    }

    return i;
}
