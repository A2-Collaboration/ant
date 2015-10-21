#include "SetupRegistry.h"

#include "Setup.h"

#include "base/Logger.h"
#include "base/std_ext/string.h"

#include <stdexcept>

using namespace std;
using namespace ant::expconfig;

SetupRegistry& SetupRegistry::get_instance()
{
    static SetupRegistry instance;
    return instance;
}

shared_ptr<Setup> SetupRegistry::GetSetup(const string& name)
{
    auto& setups = get_instance().setups;
    auto it_setup = setups.lower_bound(name);
    if(it_setup == setups.end() || setups.key_comp()(name, it_setup->first)) {
        // setup not created yet
        auto& setup_creators = get_instance().setup_creators;
        auto it_setupcreator = setup_creators.find(name);
        if(it_setupcreator == setup_creators.end())
            return nullptr;
        // found creator
        auto setup = it_setupcreator->second(name, get_instance().options);
        if(setup->GetName() != name)
            throw std::runtime_error(std_ext::formatter()
                                     << "Setup name " << name << " does not match GetName() " << setup->GetName());
        it_setup = setups.emplace_hint(it_setup, name, setup);
    }
    return it_setup->second;
}

void SetupRegistry::RegisterSetup(Creator creator, string name)
{
    setup_creators[name] = creator;
}

SetupRegistry::SetupRegistry() : options(make_shared<const OptionsList>())
{

}

SetupRegistry::~SetupRegistry()
{

}

void SetupRegistry::Destroy()
{
    get_instance().setups.clear();
}

list<string> SetupRegistry::GetNames()
{
    list<string> list;
    for(const auto& entry : get_instance().setup_creators) {
        list.emplace_back(entry.first);
    }
    return list;
}

void SetupRegistry::SetSetupOptions(SetupOptPtr opt)
{
    if(!get_instance().setups.empty())
        throw runtime_error("Set SetupOptions before any setups have been created");
    get_instance().options = opt;
}

SetupRegistration::SetupRegistration(SetupRegistry::Creator creator, string name)
{
    SetupRegistry::get_instance().RegisterSetup(creator, name);
}
