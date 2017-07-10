#pragma once

#include "Physics.h"

#include <memory>
#include <queue>

class TTree;

namespace ant {

namespace analysis {

class SlowControlManager;

namespace slowcontrol {
struct event_t;
}

namespace utils {
class ParticleID;
}

namespace input {
struct event_t;
class DataReader;
}

class PhysicsManager {
protected:
    using physics_list_t = std::list< std::unique_ptr<Physics> >;

    physics_list_t physics;

    using readers_t = std::list< std::unique_ptr<input::DataReader> >;
    readers_t amenders;
    std::unique_ptr<input::DataReader> source;

    void InitReaders(readers_t readers_);
    bool TryReadEvent(input::event_t& event);

    std::unique_ptr<SlowControlManager> slowcontrol_mgr;

    virtual void ProcessEvent(input::event_t& event, physics::manager_t& manager);
    virtual void SaveEvent(input::event_t event, const physics::manager_t& manager);

    struct interrupt_t {
        interrupt_t(volatile bool* interrupt_) :
            interrupt(interrupt_) {}
        explicit operator bool() {
            if(interrupt != nullptr)
                return *interrupt;
            return false;
        }
    private:
        volatile bool* interrupt = nullptr;
    };

    interrupt_t interrupt;

    interval<TID> processedTIDrange;

    // for output of TEvents to TTree
    TTree*  treeEvents;
    TEvent* treeEventPtr;

public:

    PhysicsManager(volatile bool* interrupt_ = nullptr);
    virtual ~PhysicsManager();

    template <typename T, typename ... args_t>
    void AddPhysics(args_t&&... args) {
       AddPhysics(
              move(std_ext::make_unique<T>(
                std::forward<args_t>(args)...
                ))
              );
    }

    void AddPhysics(std::unique_ptr<Physics> pc) {

        if(pc==nullptr)
            return;
        physics.emplace_back(std::move(pc));
    }

    const interval<TID>& GetProcessedTIDRange() const { return processedTIDrange; }

    void ReadFrom(std::list<std::unique_ptr<input::DataReader> > readers_,
                  long long maxevents
                  );

    virtual void ShowResults();

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

};

}} // namespace ant::analysis
