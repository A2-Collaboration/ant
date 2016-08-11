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
    p->Reset();
    processors.emplace_back(p);
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
            << processors.size() << " processors";
}

bool SlowControlManager::processor_t::IsComplete() {
    if(Type == type_t::Unknown)
        return false;
    if(Type == type_t::Forward)
        return true;
    // Type == Backward
    return !TIDs.empty();
}

bool SlowControlManager::ProcessEvent(TEvent event)
{
    // process the reconstructed event (if any)

    physics::manager_t manager;
    bool wants_skip = false;
    bool all_complete = true;

    for(auto& p : processors) {

        TEventData& reconstructed = event.Reconstructed();

        const auto result = p.Processor->ProcessEventData(reconstructed, manager);

        if(result == slowcontrol::Processor::return_t::Complete) {
            // forward processors become complete first before
            // they request processing events, as opposed to backward processors,
            // which buffer at least one event before becoming complete

            if(p.Type == processor_t::type_t::Unknown)
                p.Type = processor_t::type_t::Forward;
            else
                p.TIDs.push(reconstructed.ID);
        }
        else if(result == slowcontrol::Processor::return_t::Skip) {
            wants_skip = true;
        }
        else if(result == slowcontrol::Processor::return_t::Process) {
            if(p.Type == processor_t::type_t::Backward)
                throw std::runtime_error("Changing types of processor is not implemented");
            p.Type = processor_t::type_t::Forward;
        }
        else if(result == slowcontrol::Processor::return_t::Buffer) {
            if(p.Type == processor_t::type_t::Forward)
                throw std::runtime_error("Changing types of processor is not implemented");
            p.Type = processor_t::type_t::Backward;
        }

        all_complete &= p.IsComplete();
    }

    // SavedForSlowControls might already be true from previous filter runs
    // so don't reset it (best we can do here, filtering and slowcontrol stuff is tricky)
    event.SavedForSlowControls |= manager.saveEvent;

    if(!wants_skip || event.SavedForSlowControls) {
        // a skipped event could still be saved in order to trigger
        // slow control processsors (see for example AcquScalerProcessor),
        // but should NOT be processed by physics classes. Mark the event accordingly in eventbuffer
        eventbuffer.emplace(wants_skip, std::move(event));
    }

    return all_complete;
}

event_t SlowControlManager::PopEvent() {

    if(eventbuffer.empty())
        return {};

    auto& front = eventbuffer.front();
    if(front.Event.HasReconstructed()) {
        // check first if all processors are still complete
        for(auto& p : processors) {
            if(!p.IsComplete()) {
                return {};
            }
        }

        const auto& id = front.Event.Reconstructed().ID;
        for(auto& p : processors) {

            if(p.Type == processor_t::type_t::Backward) {
                if(p.TIDs.front() == id) {
                    p.TIDs.pop();
                    front.PopAfter.push_back(p.Processor);
                }
            }
            else if(p.Type == processor_t::type_t::Forward) {
                if(!p.TIDs.empty() && p.TIDs.front() == id) {
                    p.TIDs.pop();
                    p.Processor->PopQueue();
                }
            }
        }
    }

    auto event = std::move(eventbuffer.front());
    eventbuffer.pop();
    return event;
}


