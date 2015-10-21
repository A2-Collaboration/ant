#include "Physics.h"

#include "base/Logger.h"
#include <stdexcept>
#include "base/std_ext/string.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;

string Physics::GetOption(const string& key) const
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

PhysicsRegistry& PhysicsRegistry::get_instance()
{
    static PhysicsRegistry instance;
    return instance;
}

std::unique_ptr<Physics> PhysicsRegistry::Create(const string& name, PhysOptPtr opts)
{
    auto physics = PhysicsRegistry::get_instance().physics_creators.at(name)(name, opts);
    if(physics->GetName() != name)
        throw std::runtime_error(std_ext::formatter()
                                 << "Physics name " << name << " does not match GetName() " << physics->GetName());
    return physics;
}

std::vector<string> PhysicsRegistry::GetList()
{
    std::vector<std::string> list;
    for(const auto& entry : get_instance().physics_creators) {
        list.emplace_back(entry.first);
    }
    return list;
}


PhysicsRegistration::PhysicsRegistration(physics_creator c, const string& name)
{
    PhysicsRegistry::get_instance().RegisterPhysics(c,name);
}




