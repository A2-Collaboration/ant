#pragma once

#include "Physics.h"
#include "tree/TSlowControl.h"
#include "analysis/data/Slowcontrol.h"
#include "analysis/physics/SlowcontrolManager.h"


namespace ant {

struct TAntHeader;

namespace analysis {

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
    readers_t readers;
    std::unique_ptr<input::DataReader> source;

    bool InitReaders(readers_t readers_);
    bool TryReadEvent(std::unique_ptr<data::Event>& event);

    slowontrol::Manager slowcontrol_mgr;
    data::Slowcontrol slowcontrol_data;

    std::queue< std::unique_ptr<data::Event> > eventbuffer;

    long long nEventsProcessed = 0;

    void ProcessEventBuffer(long long maxevents);
    void ProcessEvent(std::unique_ptr<data::Event> event);

    bool progressUpdates = true;



    struct running_t {
        running_t(volatile bool* running_) :
            running(running_) {}
        explicit operator bool() {
            if(running != nullptr)
                return *running;
            return true;
        }
    private:
        volatile bool* running = nullptr;
    };
    running_t running;
    TID firstID;
    TID lastID;

public:

    PhysicsManager(volatile bool* running_ = nullptr);
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

        pc->Initialize(slowcontrol_data);
        physics.emplace_back(std::move(pc));
    }

    void SetParticleID(std::unique_ptr<utils::ParticleID> pid);

    void SetAntHeader(TAntHeader& header);


    void ReadFrom(std::list<std::unique_ptr<input::DataReader> > readers_,
                  long long maxevents
                  );

    void ShowResults();

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

    void EnableProgressUpdates(bool updates=false);
};

}} // namespace ant::analysis
