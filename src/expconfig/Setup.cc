#include "Setup.h"

#include "base/Logger.h"

#include <string>
#include <vector>
#include <stdexcept>

using namespace ant::expconfig;

void SetupRegistry::init_setups()
{
    if(setups.size() == setup_creators.size())
        return;
    for(const auto& creator : setup_creators) {
        const auto& setup = creator();
        VLOG(4) << "Adding config " << setup->GetName();
        setups.emplace_back(move(setup));
    }
}

SetupRegistry& SetupRegistry::get()
{
    static SetupRegistry instance;
    return instance;
}

void SetupRegistry::add(setup_creator creator)
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

SetupRegistration::SetupRegistration(setup_creator creator)
{
    SetupRegistry::get().add(creator);
}
