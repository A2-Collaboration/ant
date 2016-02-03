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
    friend class ant::analysis::PhysicsManager;
    friend class ant::analysis::SlowControlManager;
    bool saveEvent = false;
    bool keepReadHits = false;

};

}
}
}
