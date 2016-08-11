#pragma once

#include "SlowControlProcessors.h"
#include "tree/TEvent.h"
#include <list>
#include <queue>
#include <memory>
#include <functional>

namespace ant {

namespace analysis {

namespace slowcontrol {

struct event_t {
    // WantsSkip indicates event which should not be processed by physics class,
    // but could still be saved in treeEvents by PhysicsManager
    bool   WantsSkip = false;
    TEvent Event;
    event_t() {}
    event_t(bool wantsSkip, TEvent event) :
        WantsSkip(wantsSkip), Event(std::move(event))
    {}
    event_t& operator=(event_t&&) = default;
    event_t(event_t&&) = default;

    // makes "while(auto e = scm.PopEvent()) {}" loops possible
    explicit operator bool() const {
        return static_cast<bool>(Event);
    }

    ~event_t() {
        for(auto& action : DeferredActions)
            action();
    }

protected:
    // PopAfter contains processors which should change after the
    // event is destroyed. This should only be accessed by SlowControlManager.
    friend class ant::analysis::SlowControlManager;
    using action_t = std::function<void(void)>;
    std::list<action_t> DeferredActions;
};

} // namespace slowcontrol

class SlowControlManager {
public:


protected:

    std::queue<slowcontrol::event_t> eventbuffer;

    using ProcessorPtr = std::shared_ptr<slowcontrol::Processor>;

    struct processor_t {
        processor_t(ProcessorPtr proc) : Processor(proc) {}
        ProcessorPtr    Processor;
        std::list<TID>  CompletionPoints;

        enum class type_t {
            Unknown, Backward, Forward
        };
        type_t Type = type_t::Unknown;

        bool IsComplete();
    };

    std::vector<processor_t> processors;

    void AddProcessor(ProcessorPtr p);

public:
    SlowControlManager();

    bool ProcessEvent(TEvent event);

    slowcontrol::event_t PopEvent();

    size_t BufferSize() const { return eventbuffer.size(); }

};

}
}
