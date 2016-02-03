#pragma once

#include "SlowControlProcessors.h"
#include "tree/TEvent.h"
#include <map>
#include <queue>
#include <memory>

namespace ant {

namespace analysis {

class SlowControlManager {
protected:

    std::queue<TEventPtr> eventbuffer;

    using buffer_t = std::queue<TID>;
    std::map<std::shared_ptr<slowcontrol::Processor>, buffer_t> slowcontrol;
    bool all_complete = false;

    void AddProcessor(std::shared_ptr<slowcontrol::Processor> p);

    static TID min(const TID& a, const TID& b);

    TID PopBuffers(const TID& id);

    TID changepoint;

public:
    SlowControlManager();

    void ProcessEvent(TEventPtr event);

    TEventPtr PopEvent();
    bool hasEvents() const;

};

}
}
