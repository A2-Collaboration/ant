#include "SlowControlManager.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "SlowControlVariables.h"

#include "base/Logger.h"

#include <stdexcept>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::slowcontrol;


void SlowControlManager::AddProcessor(ProcessorPtr p)
{
    slowcontrol.insert({p, buffer_t()});
}

SlowControlManager::SlowControlManager()
{
    unsigned nRegistered = 0;
    for(VariablePtr var : Variables::All) {
        if(var->GetProcessors().empty())
            continue;

        nRegistered++;
        for(auto& p : var->GetProcessors())
            AddProcessor(p);
    }

    LOG_IF(nRegistered>0, INFO)
            << "Have " << nRegistered
            << " registered slowcontrol variables using "
            << slowcontrol.size() << " processors";
}

bool SlowControlManager::ProcessEvent(TEvent event)
{

    // at a changepoint, pop the slow control
    if(!changepoint.IsInvalid()) {
        for(auto& sl : slowcontrol) {
            const ProcessorPtr& p = sl.first;
            buffer_t& tid_buffer = sl.second;

            assert(!tid_buffer.empty());

            if(tid_buffer.front() == changepoint) {
                tid_buffer.pop();
                p->PopQueue();
            }
        }
        changepoint = {};
    }

    if(!event.HasReconstructed())
        return true;

    // always process the event
    {
        TEventData& reconstructed = event.Reconstructed();

        physics::manager_t manager;
        bool wants_skip = false;
        all_complete = true;

        for(auto& sl : slowcontrol) {
            const ProcessorPtr& p = sl.first;
            buffer_t& tid_buffer = sl.second;

            const auto result = p->ProcessEventData(reconstructed, manager);

            if(result == slowcontrol::Processor::return_t::Complete) {
                tid_buffer.push(reconstructed.ID);
            }
            else if(result == slowcontrol::Processor::return_t::Skip) {
                wants_skip = true;
            }

            all_complete &= !tid_buffer.empty();
        }

        // SavedForSlowControls might already be true from previous filter runs
        event.SavedForSlowControls |= manager.saveEvent;

        if(!wants_skip || event.SavedForSlowControls) {
            // a skipped event could still be saved in order to trigger
            // slow control processsors (see for example AcquScalerProcessor),
            // but should NOT be processed by physics classes. Mark the event accordingly in eventbuffer
            eventbuffer.emplace(wants_skip, std::move(event));
        }

    }

    // update the changepoint
    if(all_complete) {
        for(auto& sl : slowcontrol) {
            const ProcessorPtr& p = sl.first;
            buffer_t& tid_buffer = sl.second;

            assert(!tid_buffer.empty());

            changepoint = changepoint.IsInvalid() ?
                              tid_buffer.front() :
                              std::min(changepoint, tid_buffer.front());
        }
    }

    return all_complete;
}

event_t SlowControlManager::PopEvent() {

    if(!all_complete || eventbuffer.empty())
        return {};

    const auto& id = eventbuffer.front().Event.Reconstructed().ID;

    // at a changepoint, we still analyse the event, but stop afterwards
    if(id == changepoint)
        all_complete = false;

//    {
//        TID min_tid;
//        for(auto& sl : slowcontrol) {
//            const ProcessorPtr& p = sl.first;
//            buffer_t& tid_buffer = sl.second;

//            if(tid_buffer.front() == id) {
//                tid_buffer.pop();
//                p->PopQueue();

//            }
//            if(tid_buffer.empty())
//                all_complete = false;
//            else
//                min_tid = min_tid.IsInvalid() ?
//                              tid_buffer.front() :
//                              std::min(min_tid, tid_buffer.front());
//        }
//        changepoint = min_tid;
//    }

    auto event = std::move(eventbuffer.front());
    eventbuffer.pop();
    return event;
}
