#pragma once

namespace ant {
namespace analysis {

struct PhysicsManager;

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
    manager_t() { Reset(); }
    void Reset() { saveEvent = false; keepReadHits = false; }
};

}
}
}