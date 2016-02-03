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
    bool      Save;  // indicates that slowcontrol processors found this interesting
    TEventPtr Event;
    event_t() {}
    event_t(bool save, TEventPtr event) :
        Save(save), Event(move(event))
    {}
    explicit operator bool() const {
        return Event != nullptr;
    }
};

} // namespace slowcontrol

class SlowControlManager {
public:


protected:

    std::queue<slowcontrol::event_t> eventbuffer;

    using buffer_t = std::queue<TID>;
    std::map<std::shared_ptr<slowcontrol::Processor>, buffer_t> slowcontrol;
    bool all_complete = false;

    void AddProcessor(std::shared_ptr<slowcontrol::Processor> p);

    static TID min(const TID& a, const TID& b);

    TID PopBuffers(const TID& id);

    TID changepoint;

public:
    SlowControlManager();

    bool ProcessEvent(TEventPtr event);

    slowcontrol::event_t PopEvent();

    size_t BufferSize() const { return eventbuffer.size(); }

};

}
}
