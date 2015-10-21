#pragma once

#include <map>
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
    using setup_creators_t = std::map<std::string, Creator>;
    using setups_t = std::map<std::string, std::shared_ptr<Setup> >;
    setup_creators_t setup_creators;
    setups_t setups;
    void RegisterSetup(Creator, std::string);
    SetupRegistry() {}
    static SetupRegistry& get_instance();
public:
    static std::shared_ptr<Setup> GetSetup(const std::string& name);
    static std::list<std::string> GetNames();
    static void Destroy();
};

/**
 * @brief The SetupRegistration class is instantiated for each setup
 */
class SetupRegistration
{
public:
    SetupRegistration(SetupRegistry::Creator, std::string);
};

#define AUTO_REGISTER_SETUP(setup) \
    SetupRegistration _setup_registration_ ## setup(std::make_shared<setup>, #setup);

}} // namespace ant::expconfig
