#include "Physics.h"

#include "base/Logger.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;

ant::analysis::OptionsList PhysicsOptions;

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



OptionsList::OptionsList(std::shared_ptr<const OptionsList> Parent):
    parent(Parent)
{
}

void OptionsList::SetOption(const string& str)
{
    const std::string delimiter = "=";
    const auto delimiter_pos = str.find(delimiter);
    if( delimiter_pos != str.npos) {
        const std::string key = str.substr(0, delimiter_pos);
        const std::string val = str.substr(delimiter_pos + delimiter.length(), str.npos);
        options.insert({key,val});
    } else {
        LOG(WARNING) << "Can't parse option string \"" << str << "\"";
    }
}

void OptionsList::SetOptions(const string& str)
{
    string::size_type p = 0;
    string::size_type np = 0;

    do {
        np = str.find(",", p);
        const auto o = str.substr(p,np);
        SetOption(o);
        p = np+1;
    } while(np != str.npos);
}

string OptionsList::GetOption(const string& key) const
{
    const auto entry = options.find(key);

    if(entry == options.end()) {
        if(parent) {
            return parent->GetOption(key);
        }
        return "";
    }

    return entry->second;

}
