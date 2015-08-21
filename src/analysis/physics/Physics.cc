#include "Physics.h"

#include "base/Logger.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;

ant::analysis::OptionsList PhysicsOptions;

const string Physics::GetOption(const string& key) const
{

    return PhysicsOptions.GetOption(key);
}

Physics::Physics(const string &name):
    name_(name),
    HistFac(name)
{}

PhysicsRegistry& PhysicsRegistry::get()
{
    static PhysicsRegistry instance;
    return instance;
}

std::unique_ptr<Physics> PhysicsRegistry::Create(const string& name)
{
    return PhysicsRegistry::get().physics_creators.at(name)();
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



void OptionsList::SetOption(const string& str)
{
    const std::string delimiter = "=";
    const auto delimiter_pos = str.find(delimiter);
    if( delimiter_pos != str.npos) {
        const std::string key = str.substr(0, delimiter_pos);
        const std::string val = str.substr(delimiter_pos + delimiter.length(), str.npos);
        options.insert({key,val});
        VLOG(9) << "Added Physics option " << key << " = " << val;
    } else {
        LOG(WARNING) << "Can't parse option string \"" << str << "\"";
    }
}

string OptionsList::GetOption(const string& key) const
{
    const auto entry = options.find(key);

    if(entry == options.end()) {
        return "";
    }

    return entry->second;

}
