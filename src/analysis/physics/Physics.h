#pragma once

#include "analysis/plot/HistogramFactories.h"
#include "analysis/input/slowcontrol/SlowControl.h"

#include "tree/TEvent.h"

#include "base/OptionsList.h"

#include <list>
#include <string>
#include <map>
#include <memory>
#include <functional>

class TFile;
class TDirectory;

namespace ant {

namespace analysis {

using PhysOptPtr = std::shared_ptr<const OptionsList>;

class Physics {
private:
    std::string name_;

protected:
    SmartHistFactory HistFac;
    const PhysOptPtr Options;

public:
    Physics(const std::string& name, PhysOptPtr opts);
    virtual ~Physics() {}

    struct manager_t {
        friend class PhysicsManager;

        void SaveEvent() {
            saveEvent = true;
        }
    private:
        manager_t() : saveEvent(false) {}
        bool saveEvent;
    };

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) =0;
    virtual void Finish() {}
    virtual void ShowResult() {}
    std::string GetName() const { return name_; }

    virtual void Initialize(input::SlowControl& slowcontrol);

    Physics(const Physics&) = delete;
    Physics& operator=(const Physics&) = delete;
};



using physics_creator = std::function<std::unique_ptr<Physics>(const std::string&, PhysOptPtr)>;

class PhysicsRegistry
{
    friend class PhysicsRegistration;

private:
    using physics_creators_t = std::map<std::string, physics_creator>;
    physics_creators_t physics_creators;
    static PhysicsRegistry& get_instance();

    void RegisterPhysics(physics_creator c, const std::string& name) {
        physics_creators[name] = c;
    }
public:

    static std::unique_ptr<Physics> Create(const std::string& name, PhysOptPtr opts = std::make_shared<OptionsList>());

    static std::vector<std::string> GetList();

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

};

class PhysicsRegistration
{
public:
    PhysicsRegistration(physics_creator c, const std::string& name);
};

template<class T>
std::unique_ptr<Physics> physics_factory(const std::string& name, PhysOptPtr opts)
{
    return std_ext::make_unique<T>(name, opts);
}

#define AUTO_REGISTER_PHYSICS(physics) \
    ant::analysis::PhysicsRegistration _physics_registration_ ## physics(ant::analysis::physics_factory<physics>, #physics);

}} // nammespace ant::analysis
