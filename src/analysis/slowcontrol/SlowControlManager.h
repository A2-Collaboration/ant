#pragma once

#include "analysis/slowcontrol/SlowControl.h"

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

    void SetRequiredKeys(const std::list<TSlowControl::Key> keys);

    void ProcessSlowControls(TEvent& event);

    /**
     * @brief check if at least one slow control variable has been requrested
     * @return
     */
    bool hasRequests() const { return !slowcontrol.empty(); }

    bool isComplete() const;

    TID FindMinimalTID() const;

    TID UpdateSlowcontrolData(slowcontrol::SlowControl& slc);

};

}
}
