#include "Physics.h"

#include "base/Logger.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;

const string Physics::GetOption(const string& key) const
{
    if(options) {
        return options->GetOption(key);
    }

    return "";
}

Physics::Physics(const string &name, PhysOptPtr opts):
    name_(name),
    options(opts),
    HistFac(name)
{}

void Physics::Initialize(data::Slowcontrol&)
{
}

PhysicsRegistry& PhysicsRegistry::get()
{
    static PhysicsRegistry instance;
    return instance;
}

std::unique_ptr<Physics> PhysicsRegistry::Create(const string& name, PhysOptPtr opts)
{
    return PhysicsRegistry::get().physics_creators.at(name)(opts);
}

std::vector<string> PhysicsRegistry::GetList() const
{
    std::vector<std::string> list(physics_creators.size());
    list.resize(0);

    for(const auto& entry : physics_creators) {
        list.emplace_back(entry.first);
    }

    return list;
}

void PhysicsRegistry::PrintRegistry()
{
    for(auto& entry : get().physics_creators) {
        LOG(INFO) << entry.first;
    }
}


PhysicsRegistration::PhysicsRegistration(physics_creator c, const string& name)
{
    PhysicsRegistry::get().RegisterPhysics(c,name);
}




