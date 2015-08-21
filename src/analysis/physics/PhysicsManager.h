#pragma once

#include "Physics.h"

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

public:
    std::shared_ptr<utils::ParticleID> particleID;

    PhysicsManager();
    template <typename T, typename ... args_t>
    void AddPhysics(args_t&&... args) {
        physics.push_back(
              std_ext::make_unique<T>(
                std::forward<args_t>(args)...
                )
              );
    }

    void AddPhysics(std::unique_ptr<Physics> pc) {
        if(pc==nullptr)
            return;
        physics.emplace_back(std::move(pc));
    }

    void ReadFrom(std::list< std::unique_ptr<input::DataReader> > readers,
                  long long maxevents,
                  bool& running,
                  TAntHeader* header
                  );


    void ProcessEvent(data::Event& event);
    void ShowResults();

};

}} // namespace ant::analysis
