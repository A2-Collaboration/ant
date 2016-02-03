#pragma once

#include "physics/manager_t.h"
#include "tree/TSlowControl.h"

#include <map>
#include <queue>
#include <memory>
#include <list>

namespace ant {

struct TEvent;

namespace analysis {

class SlowControlManager {
protected:

    using buffer_t = std::queue<std::pair<TID,  TSlowControl>>;
    std::map<TSlowControl::Key, buffer_t> slowcontrol;

    TID minimal_in_buffer;

    static TID min(const TID& a, const TID& b);

public:
    SlowControlManager();

    void ProcessEvent(const TEvent& event, physics::manager_t& manager);

    /**
     * @brief check if at least one slow control variable has been requrested
     * @return
     */
    bool hasRequests() const { return !slowcontrol.empty(); }

    bool isComplete() const;

    TID GetRunUntil() const;

};

}
}
