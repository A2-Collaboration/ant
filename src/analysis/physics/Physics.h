#pragma once

#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"

#include "base/std_ext.h"

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


struct OptionsList {
    std::map<std::string, std::string> options;

    void SetOption(const std::string& str);

    std::string GetOption(const std::string& key) const;

};

using PhysOptPtr = std::shared_ptr<const OptionsList>;

class Physics {
private:
    std::string name_;
    PhysOptPtr options;

protected:
    SmartHistFactory HistFac;

    const std::string GetOption(const std::string& key) const;

public:
    Physics(const std::string& name, PhysOptPtr opts=nullptr);
    virtual ~Physics() {}
    virtual void ProcessEvent(const data::Event& event) =0;
    virtual void Finish() =0;
    virtual void ShowResult() =0;
    std::string GetName() { return name_; }
};

template<class T>
std::unique_ptr<Physics> physics_factory(PhysOptPtr opts)
{
    return std::move(std_ext::make_unique<T>(opts));
}

using physics_creator = std::function<std::unique_ptr<Physics>(PhysOptPtr)>;

class PhysicsRegistry
{
private:
    using physics_creators_t = std::map<std::string, physics_creator>;
    physics_creators_t physics_creators;

public:
    static PhysicsRegistry& get();

    static std::unique_ptr<Physics> Create(const std::string& name, PhysOptPtr opts);

    std::vector<std::string> GetList() const;

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
    ant::analysis::PhysicsRegistration _physics_registration_ ## physics(ant::analysis::physics_factory<physics>,name);

}} // nammespace ant::analysis
