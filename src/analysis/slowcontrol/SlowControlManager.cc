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
    p->Init();
    processors.emplace_back(p);
}

SlowControlManager::SlowControlManager(const input::reader_flags_t& reader_flags)
{
    unsigned nRegistered = 0;
    for(VariablePtr var : Variables::All) {
        if(!var->requested)
            continue;

        nRegistered++;

        var->Init(reader_flags);
        for(auto& p : var->GetNeededProcessors())
            AddProcessor(p);
    }

    LOG_IF(nRegistered>0, INFO)
            << "Have " << nRegistered
            << " registered slowcontrol variables using "
            << processors.size() << " processors";
}

bool SlowControlManager::processor_t::IsComplete() const {
    if(Type == type_t::Unknown)
        return false;
    return !CompletionPoints.empty();
}

bool SlowControlManager::ProcessEvent(input::event_t event)
{
    // process the reconstructed event (if any)

    physics::manager_t manager;
    bool wants_skip = false;
    bool all_complete = true;

    for(auto& p : processors) {

        TEventData& reconstructed = event.Reconstructed();

        const auto result = p.Processor->ProcessEventData(reconstructed, manager);

        if(result == slowcontrol::Processor::return_t::Complete) {
            p.CompletionPoints.push_back(reconstructed.ID);
        }
        else if(result == slowcontrol::Processor::return_t::Skip) {
            wants_skip = true;
        }
        else if(result == slowcontrol::Processor::return_t::Process) {
            if(p.Type == processor_t::type_t::Backward)
                throw std::runtime_error("Backward+Forward processor is not supported");
            p.Type = processor_t::type_t::Forward;
        }
        else if(result == slowcontrol::Processor::return_t::Buffer) {
            if(p.Type == processor_t::type_t::Forward)
                throw std::runtime_error("Backward+Forward processor is not supported");
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

slowcontrol::event_t SlowControlManager::PopEvent() {

    if(eventbuffer.empty())
        return {};

    auto& front = eventbuffer.front();
    if(front.Event.HasReconstructed()) {


        // check first if all processors are still complete
        // otherwise go back to filling
        for(auto& p : processors)
            if(!p.IsComplete())
                return {};

        // check for potential pop/change points
        const auto& id = front.Event.Reconstructed().ID;
        for(auto& p : processors) {
            auto proc = p.Processor;

            if(p.Type == processor_t::type_t::Backward) {
                if(p.CompletionPoints.front() == id) {
                    p.CompletionPoints.pop_front();
                    front.DeferredActions.emplace_back(
                                [proc] () {
                        proc->PopQueue();
                        proc->SetHasChanged(true);
                    });

                }
            }
            else if(p.Type == processor_t::type_t::Forward) {
                // handle changing the processor's value
                if(p.CompletionPoints.size()>1) {
                    auto second_to_front = std::next(p.CompletionPoints.begin());
                    if(*second_to_front == id) {
                        p.CompletionPoints.pop_front();
                        proc->PopQueue();
                        proc->SetHasChanged(true);
                    }
                }
            }

            // upkeep HasChanged until non-skipped event is popped
            // this works for forward/backward
            if(proc->HasChanged() && !front.WantsSkip) {
                front.DeferredActions.emplace_back(
                            [proc] () {
                    proc->SetHasChanged(false);
                });
            }
        }
    }

    auto event = std::move(eventbuffer.front());
    eventbuffer.pop();
    return event;
}


