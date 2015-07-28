#include "ExpConfig.h"

#include "Setup.h"
#include "tree/THeaderInfo.h"
#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"
#include "base/Logger.h"

#include <type_traits>
#include <list>
#include <iostream>

using namespace std;

namespace ant { // template implementations need explicit namespace

std::string ExpConfig::ManualSetupName = ""; // default empty

template<typename T>
shared_ptr<T> ExpConfig::Get_(const THeaderInfo& header) {
    static_assert(is_base_of<Base, T>::value, "T must be a base of ExpConfig::Base");

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

    shared_ptr<Base> config = nullptr;

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
        // automatic search mode in all registerd configs

        // make a copy of the list of registered configs
        // they only need to implement the Matches(THeaderInfo&) methodd
        auto& registry = expconfig::SetupRegistry::get();
        std::list< std::shared_ptr<Base> > modules(registry.begin(), registry.end());

        // remove the config if the config says it does not match
        modules.remove_if([&header] (const shared_ptr<Base>& m) {
            return !m->Matches(header);
        });

        // check if something reasonable is left
        if(modules.empty()) {
            throw ExpConfig::Exception(std_ext::formatter()
                                       << "No config found for header "
                                   << header);
        }
        if(modules.size()>1) {
            throw ExpConfig::Exception(std_ext::formatter()
                                       << "More than one config found for header "
                                       << header);
        }
        config = modules.back();
    }

    // try to cast the found config instance
    const auto& ptr = dynamic_pointer_cast<T, Base>(config);

    if(ptr==nullptr) {
        throw Exception(
                    std_ext::formatter()
                    << "Found config does not fit to requested type "
                    << std_ext::getTypeAsString<T>());
    }

    // hand over the unique ptr
    return ptr;
}


shared_ptr<ExpConfig::Setup> ExpConfig::Setup::Get(const THeaderInfo& header)
{
    return Get_<ExpConfig::Setup>(header);
}

list< shared_ptr<ExpConfig::Setup> > ExpConfig::Setup::GetAll() {
    // make a copy of the list of registered configs
    auto& registry = expconfig::SetupRegistry::get();
    std::list< std::shared_ptr<ExpConfig::Setup> > modules(registry.begin(), registry.end());
    return modules;
}

shared_ptr<ExpConfig::Setup> ExpConfig::Setup::Get(const std::string& name)
{
    for(const auto& module : GetAll()) {
        if(module->GetName() == name)
            return module;
    }
    // nothing found
    return nullptr;
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




