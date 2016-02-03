#pragma once

namespace ant {
namespace analysis {

class PhysicsManager;
class SlowControlManager;

namespace physics {

struct manager_t {
    void SaveEvent() {
        saveEvent = true;
    }
    void KeepDetectorReadHits() {
        keepReadHits = true;
    }
private:
    bool saveEvent;
    bool keepReadHits;
    friend class ant::analysis::PhysicsManager;
    friend class ant::analysis::SlowControlManager;
    manager_t() { Reset(); }
    void Reset() { saveEvent = false; keepReadHits = false; }
};

}
}
}
