#include "SetupRegistry.h"

#include "Setup.h"

#include "base/Logger.h"

using namespace ant::expconfig;

void SetupRegistry::init_setups()
{
    if(setups.size() == setup_creators.size())
        return;
    for(const auto& creator : setup_creators) {
        const auto& setup = creator();
        VLOG(9) << "Adding setup config " << setup->GetName() << " to registry";
        setups.emplace_back(move(setup));
    }
}

SetupRegistry& SetupRegistry::get()
{
    static SetupRegistry instance;
    return instance;
}

void SetupRegistry::add(Creator creator)
{
    setup_creators.push_back(creator);
}

SetupRegistry::setups_t::iterator SetupRegistry::begin()
{
    init_setups();
    return setups.begin();
}

SetupRegistry::setups_t::iterator SetupRegistry::end()
{
    init_setups();
    return setups.end();
}

void SetupRegistry::destroy()
{
    setups.clear();
}

SetupRegistration::SetupRegistration(SetupRegistry::Creator creator)
{
    SetupRegistry::get().add(creator);
}
