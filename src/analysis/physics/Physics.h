#pragma once

#include "analysis/data/Event.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/input/DataReader.h"

#include "base/std_ext.h"

#include <list>
#include <string>
#include <memory>
#include <functional>

class TFile;
class TDirectory;

namespace ant {

class Physics {
protected:
    SmartHistFactory HistFac;
public:
    Physics(const std::string& name);
    virtual ~Physics() {}
    virtual void ProcessEvent(const ant::Event& event) =0;
    virtual void Finish() =0;
    virtual void ShowResult() =0;
};

class DebugPhysics: public Physics {
public:
    DebugPhysics(const std::string& name="DebugPhysics"): Physics(name) {}
    virtual ~DebugPhysics() {}

    virtual void ProcessEvent(const ant::Event& event);
    virtual void Finish();
    virtual void ShowResult();
};

class PhysicsManager {
protected:
    using physics_list_t = std::list< std::unique_ptr<ant::Physics> >;

    physics_list_t physics;

public:
    PhysicsManager();
    template <typename T, typename ... args_t>
    void AddPhysics(args_t&&... args) {
        physics.push_back(
              std_ext::make_unique<T>(
                std::forward<args_t>(args)...
                )
              );
    }

    void AddPhysics(std::unique_ptr<ant::Physics> pc) {
        if(pc==nullptr)
            return;
        physics.emplace_back(std::move(pc));
    }

    void ReadFrom(std::list< std::unique_ptr<input::DataReader> > readers,
                  long long maxevents,
                  bool& running
                  );


    void ProcessEvent(const ant::Event& event);
    void ShowResults();

};

template<class T>
std::unique_ptr<Physics> physics_factory()
{
    return std::move(std_ext::make_unique<T>());
}

using physics_creator = std::function<std::unique_ptr<ant::Physics>()>;

class PhysicsRegistry
{
private:
    using physics_creators_t = std::map<std::string, physics_creator>;
    physics_creators_t physics_creators;

public:
    static PhysicsRegistry& get();

    static std::unique_ptr<ant::Physics> Create(const std::string& name);

    void RegisterPhysics(physics_creator c, const std::string& name) {
        physics_creators[name] = c;
    }

    static void PrintRegistry();

};

class PhysicsRegistration
{
public:
    PhysicsRegistration(physics_creator c, const std::string& name);
};

#define AUTO_REGISTER_PHYSICS(physics, name) \
    ant::PhysicsRegistration _physics_registration_ ## physics(ant::physics_factory<physics>,name);

}
