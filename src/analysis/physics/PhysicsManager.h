#pragma once

#include "Physics.h"

#include <memory>
#include <queue>

class TTree;

namespace ant {

struct TAntHeader;
struct TEvent;

namespace analysis {

class SlowControlManager;

namespace slowcontrol {
struct event_t;
}

namespace utils {
class ParticleID;
}

namespace input {
class DataReader;
}

class PhysicsManager {
protected:
    using physics_list_t = std::list< std::unique_ptr<Physics> >;

    physics_list_t physics;

    std::unique_ptr<utils::ParticleID> particleID;

    using readers_t = std::list< std::unique_ptr<input::DataReader> >;
    readers_t amenders;
    std::unique_ptr<input::DataReader> source;

    void InitReaders(readers_t readers_);
    bool TryReadEvent(TEvent& event);

    std::unique_ptr<SlowControlManager> slowcontrol_mgr;

    virtual void ProcessEvent(TEvent& event, physics::manager_t& manager);
    virtual void SaveEvent(slowcontrol::event_t event, const physics::manager_t& manager);

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

    TID firstID;
    TID lastID;

    // for output of TEvents to TTree
    TTree*  treeEvents;
    TEvent* treeEventPtr;

    // for progress
    double last_PercentDone;

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

    void SetParticleID(std::unique_ptr<utils::ParticleID> pid);

    void SetAntHeader(TAntHeader& header);


    void ReadFrom(std::list<std::unique_ptr<input::DataReader> > readers_,
                  long long maxevents
                  );

    virtual void ShowResults();

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

};

}} // namespace ant::analysis
