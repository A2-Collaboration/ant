#pragma once

#include "Physics.h"
#include "tree/TSlowControl.h"
#include "analysis/physics/slowcontrol/Slowcontrol.h"
#include "analysis/data/Slowcontrol.h"

#include <queue>

namespace ant {

class TAntHeader;

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

    std::map<TSlowControl::Key, std::queue< std::unique_ptr<TSlowControl> > > slowcontrol;

    std::queue< std::unique_ptr<data::Event> > eventbuffer;

    long long nEventsProcessed = 0;

    void ProcessEventBuffer(long long maxevents, bool& running, TAntHeader& header);
    void ProcessEvent(std::unique_ptr<data::Event> event);

    ant::analysis::slowcontrol::Distributor slowcontrolDistributor;

    ant::analysis::data::Slowcontrol slc;

public:

    PhysicsManager();
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

        pc->Initialize(slc);
        physics.emplace_back(std::move(pc));
    }

    void SetParticleID(std::unique_ptr<utils::ParticleID> pid);

    void ReadFrom(std::list<std::unique_ptr<input::DataReader> > readers_,
                  long long maxevents,
                  bool& running,
                  TAntHeader& header
                  );

    void ShowResults();
};

}} // namespace ant::analysis
