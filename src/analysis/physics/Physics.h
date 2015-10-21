#pragma once

#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"

#include "analysis/data/Slowcontrol.h"

#include "base/OptionsList.h"

#include <list>
#include <string>
#include <map>
#include <memory>
#include <functional>

class TFile;
class TDirectory;

namespace ant {

class TAntHeader;

namespace analysis {

using PhysOptPtr = std::shared_ptr<const OptionsList>;

class Physics {
private:
    std::string name_;
    PhysOptPtr options;

protected:
    SmartHistFactory HistFac;

    std::string GetOption(const std::string& key) const;

public:
    Physics(const std::string& name, PhysOptPtr opts=nullptr);
    virtual ~Physics() {}
    virtual void ProcessEvent(const data::Event& event) =0;
    virtual void Finish() {}
    virtual void ShowResult() =0;
    std::string GetName() const { return name_; }

    virtual void Initialize(data::Slowcontrol& slowcontrol);

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

    static std::unique_ptr<Physics> Create(const std::string& name, PhysOptPtr opts);

    static std::vector<std::string> GetList();

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
