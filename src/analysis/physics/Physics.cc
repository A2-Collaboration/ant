#include "Physics.h"

#include "base/Logger.h"
#include <stdexcept>
#include "base/std_ext/string.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;

Physics::Physics(const string &name, OptionsPtr opts):
    name_(name),
    HistFac(name)
{
    if(opts)
        HistFac.SetDirDescription(opts->Flatten());
}

PhysicsRegistry& PhysicsRegistry::get_instance()
{
    static PhysicsRegistry instance;
    return instance;
}

std::unique_ptr<Physics> PhysicsRegistry::Create(const string& name, OptionsPtr opts)
{
    auto& creators = PhysicsRegistry::get_instance().physics_creators;

    auto creator = creators.find(name);
    if(creator == creators.end())
        throw std::runtime_error("Physics class " + name + " not found");

    // this may throw an exception
    auto physics = creator->second(name, opts);

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




