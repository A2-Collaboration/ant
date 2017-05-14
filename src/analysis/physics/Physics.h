#pragma once

#include "analysis/physics/manager_t.h"

// always needed by physics classes
#include "analysis/plot/HistogramFactory.h"
#include "analysis/plot/RootDraw.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"
#include "base/OptionsList.h"
#include "base/std_ext/memory.h"

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
    HistogramFactory HistFac;

public:
    Physics(const std::string& name, OptionsPtr opts);
    virtual ~Physics() {}

    virtual void ProcessEvent(const TEvent& event, physics::manager_t& manager) =0;
    virtual void Finish() {}
    virtual void ShowResult() {}
    std::string GetName() const { return name_; }

    Physics(const Physics&) = delete;
    Physics& operator=(const Physics&) = delete;

    // provide some commonly used exception classes
    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    struct ExceptionOptionNeeded : Exception {
        using Exception::Exception;
    };
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
