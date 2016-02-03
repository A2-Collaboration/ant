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

void SlowControlManager::ProcessEvent(TEventPtr event)
{
    if(!event->Reconstructed)
        return;

    auto& reconstructed = *(event->Reconstructed);

    physics::manager_t manager;

    all_complete = true;

    for(auto& sl : slowcontrol) {

        const auto result = sl.first->ProcessEventData(reconstructed, manager);

        if(result == slowcontrol::Processor::return_t::Complete) {
            sl.second.push(reconstructed.ID);
        }

        all_complete &= !sl.second.empty();
    }


    //if manager says save-> mark event to be saved

    eventbuffer.emplace(move(event));
}

TEventPtr SlowControlManager::PopEvent() {

    if(!all_complete || eventbuffer.empty())
        return nullptr;

    auto i = move(eventbuffer.front());
    eventbuffer.pop();

    if(i->Reconstructed->ID == changepoint) {
        changepoint = PopBuffers(changepoint);
    }

    return i;
}
