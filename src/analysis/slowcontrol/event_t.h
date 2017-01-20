#pragma once

#include "tree/TEvent.h"

#include <list>
#include <functional>

namespace ant {
namespace analysis {

class SlowControlManager;

namespace slowcontrol {

struct event_t {
    // WantsSkip indicates event which should not be processed by physics class,
    // but could still be saved in treeEvents by PhysicsManager
    const bool WantsSkip = false;
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
        for(const auto& action : DeferredActions)
            action();
    }

protected:
    // PopAfter contains processors which should change after the
    // event is destroyed. This should only be accessed by SlowControlManager.
    friend class ant::analysis::SlowControlManager;
    using action_t = std::function<void(void)>;
    std::list<action_t> DeferredActions;
};

}}} // namespace ant::analysis::slowcontrol