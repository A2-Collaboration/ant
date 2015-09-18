#include "ExpConfig.h"

#include "setups/Setup.h"
#include "tree/THeaderInfo.h"
#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"
#include "base/Logger.h"
#include "base/std_ext/misc.h"

#include <type_traits>
#include <list>
#include <iostream>

using namespace std;

namespace ant { // template implementations need explicit namespace

std::string ExpConfig::ManualSetupName = ""; // default empty
std::shared_ptr<ExpConfig::Setup> ExpConfig::lastSetupFound = nullptr; // default nothing found so far

template<typename T>
shared_ptr<T> ExpConfig::Get_(const THeaderInfo& header) {

    VLOG(9) << "Searching for config of type "
            << std_ext::getTypeAsString<T>();

    // check first if some manual override of the automatic
    // search is requested
    // remember if the header already contains such an information
    if(ManualSetupName.empty() && !header.SetupName.empty()) {
        ManualSetupName = header.SetupName;
        LOG(INFO) << "Manually set setup name " << ManualSetupName
                  << " obtained from header " << header;
    }

    shared_ptr<Setup> config = nullptr;

    if(!ManualSetupName.empty()) {
        config = ExpConfig::Setup::Get(ManualSetupName);
        if(config == nullptr) {
            throw Exception(
                        std_ext::formatter()
                        << "Found no config matching name "
                        << ManualSetupName
                        );
        }
    }
    else {
        // go to automatic search mode in all registered setups

        // make a copy of the list of registered configs
        // they only need to implement the Matches(THeaderInfo&) methodd
        auto& registry = expconfig::SetupRegistry::get();
        std::list< std::shared_ptr<Setup> > modules(registry.begin(), registry.end());

        // remove the config if the config says it does not match
        modules.remove_if([&header] (const shared_ptr<Setup>& m) {
            return !m->Matches(header);
        });

        // check if something reasonable is left
        if(modules.empty()) {
            throw ExpConfig::Exception(std_ext::formatter()
                                       << "No setup found for header "
                                   << header);
        }
        if(modules.size()>1) {
            throw ExpConfig::Exception(std_ext::formatter()
                                       << "More than one setup found for header "
                                       << header);
        }
        lastSetupFound = modules.back();
        config = modules.back();
    }

    // try to cast the found config instance
    const auto& ptr = dynamic_pointer_cast<T, Setup>(config);

    if(ptr==nullptr) {
        throw Exception(
                    std_ext::formatter()
                    << "Found config does not fit to requested type "
                    << std_ext::getTypeAsString<T>());
    }

    // hand over the ptr
    return ptr;
}


shared_ptr<ExpConfig::Setup> ExpConfig::Setup::Get(const THeaderInfo& header)
{
     return Get_<ExpConfig::Setup>(header);
}

list<shared_ptr<ExpConfig::Setup>> ExpConfig::Setup::getAll()
{
    auto& registry = expconfig::SetupRegistry::get();
    std::list< std::shared_ptr<ExpConfig::Setup> > modules(registry.begin(), registry.end());
    return modules;
}

shared_ptr<ExpConfig::Setup> ExpConfig::Setup::Get(const std::string& name)
{
    for(const auto& module : getAll()) {
        if(module->GetName() == name) {
            lastSetupFound = module;
            return module;
        }
    }
    // nothing found matching the name
    return nullptr;
}

shared_ptr<ExpConfig::Setup> ExpConfig::Setup::GetLastFound()
{
    // try if we have a name
    if(lastSetupFound==nullptr && !ManualSetupName.empty()) {
        Get(ManualSetupName);
    }
    return lastSetupFound;
}

std::shared_ptr<Detector_t> ExpConfig::Setup::GetDetector(Detector_t::Type_t type)
{
    auto config = std::dynamic_pointer_cast<Reconstruct, Setup>(GetLastFound());
    if(config == nullptr)
        return nullptr;
    for(const auto& detector : config->GetDetectors()) {
        if(detector->Type == type)
            return detector;
    }
    return nullptr;
}

void ExpConfig::Setup::Cleanup()
{
    lastSetupFound = nullptr;
    expconfig::SetupRegistry::get().destroy();
}

list<string> ExpConfig::Setup::GetNames() {
    list<string> names;
    for(auto setup : getAll()) {
        names.push_back(setup->GetName());
    }
    return names;
}

shared_ptr<ExpConfig::Reconstruct> ExpConfig::Reconstruct::Get(const THeaderInfo& header)
{
    return Get_<ExpConfig::Reconstruct>(header);
}

template<>
shared_ptr< UnpackerAcquConfig >
ExpConfig::Unpacker<UnpackerAcquConfig>::Get(const THeaderInfo& header) {
    return Get_< UnpackerAcquConfig >(header);
}

template<>
shared_ptr< UnpackerA2GeantConfig >
ExpConfig::Unpacker<UnpackerA2GeantConfig>::Get(const THeaderInfo& header) {

    return Get_< UnpackerA2GeantConfig >(header);
}

} // namespace ant




