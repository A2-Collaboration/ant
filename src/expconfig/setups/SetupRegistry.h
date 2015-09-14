#pragma once

#include <list>
#include <memory>
#include <functional>

namespace ant {
namespace expconfig {

class Setup;

/**
 * @brief The SetupRegistry class semi-automatically registers Setups
 *
 * \note don't forget to use AUTO_REGISTER_SETUP
 * \note linking order is important to get this registry properly working, see CMakeLists.txt
 */
class SetupRegistry
{
friend class SetupRegistration;

private:
    using Creator = std::function<std::shared_ptr<Setup>()>;
    using setup_creators_t = std::list<Creator>;
    using setups_t = std::list< std::shared_ptr<Setup> >;
    setup_creators_t setup_creators;
    setups_t setups;
    void init_setups();
    void add(Creator);
    SetupRegistry() {}
public:
    static SetupRegistry& get();
    setups_t::iterator begin();
    setups_t::iterator end();
    void destroy();
};

/**
 * @brief The SetupRegistration class is instantiated for each setup
 */
class SetupRegistration
{
public:
    SetupRegistration(SetupRegistry::Creator);
};

#define AUTO_REGISTER_SETUP(setup) \
    SetupRegistration _setup_registration_ ## setup(std::make_shared<setup>);

}} // namespace ant::expconfig
