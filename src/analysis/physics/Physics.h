#pragma once

#include "analysis/plot/HistogramFactories.h"
#include "analysis/slowcontrol/SlowControl.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

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

class Physics {
private:
    std::string name_;

protected:
    SmartHistFactory HistFac;
    const OptionsPtr Options;

public:
    Physics(const std::string& name, OptionsPtr opts);
    virtual ~Physics() {}

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
        friend class PhysicsManager;
        manager_t() { Reset(); }
        void Reset() { saveEvent = false; keepReadHits = false; }
    };

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) =0;
    virtual void Finish() {}
    virtual void ShowResult() {}
    std::string GetName() const { return name_; }

    virtual void Initialize(slowcontrol::SlowControl& slowcontrol);

    Physics(const Physics&) = delete;
    Physics& operator=(const Physics&) = delete;
};



using physics_creator = std::function<std::unique_ptr<Physics>(const std::string&, OptionsPtr)>;

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

    static std::unique_ptr<Physics> Create(const std::string& name, OptionsPtr opts = std::make_shared<OptionsList>());

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
std::unique_ptr<Physics> physics_factory(const std::string& name, OptionsPtr opts)
{
    return std_ext::make_unique<T>(name, opts);
}

#define AUTO_REGISTER_PHYSICS(physics) \
    ant::analysis::PhysicsRegistration _physics_registration_ ## physics(ant::analysis::physics_factory<physics>, #physics);

}} // nammespace ant::analysis
