#pragma once

#include "SlowControlProcessors.h"
#include "tree/TEvent.h"
#include <map>
#include <queue>
#include <memory>

namespace ant {

namespace analysis {

namespace slowcontrol {

struct event_t {
    bool   WantsSkip;  // indicates event which should not be processed by physics
    TEvent Event;
    event_t() {}
    event_t(bool wantsSkip, TEvent event) :
        WantsSkip(wantsSkip), Event(std::move(event))
    {}
    explicit operator bool() const {
        return static_cast<bool>(Event);
    }
};

} // namespace slowcontrol

class SlowControlManager {
public:


protected:

    std::queue<slowcontrol::event_t> eventbuffer;

    using buffer_t = std::queue<TID>;
    using ProcessorPtr = std::shared_ptr<slowcontrol::Processor>;
    std::map<ProcessorPtr, buffer_t> slowcontrol;

    void AddProcessor(ProcessorPtr p);


    bool all_complete = false;
    TID changepoint;

public:
    SlowControlManager();

    bool ProcessEvent(TEvent event);

    slowcontrol::event_t PopEvent();

    size_t BufferSize() const { return eventbuffer.size(); }

};

}
}
