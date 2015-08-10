#pragma once

#include <list>
#include <memory>
#include <functional>

namespace ant {
namespace expconfig {

class Setup;

// stuff for semiauto-registering the setups
// don't forget to use AUTO_REGISTER_SETUP
// note that linking order is important to get this registry
// properly working, see CMakeLists.txt

template<class T>
std::shared_ptr<Setup> setup_factory()
{
    return std::move(std::make_shared<T>());
}

using setup_creator = std::function<std::shared_ptr<Setup>()>;

class SetupRegistry
{
private:
    using setup_creators_t = std::list<setup_creator>;
    using setups_t = std::list< std::shared_ptr<Setup> >;
    setup_creators_t setup_creators;
    setups_t setups;
    void init_setups();
public:
    static SetupRegistry& get();
    void add(setup_creator);
    setups_t::iterator begin();
    setups_t::iterator end();
    void destroy();
};

class SetupRegistration
{
public:
    SetupRegistration(setup_creator);
};

#define AUTO_REGISTER_SETUP(setup) \
    SetupRegistration _setup_registration_ ## setup(setup_factory<setup>);

}} // namespace ant::expconfig
