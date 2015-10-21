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

std::string ExpConfig::Setup::ManualName = ""; // default empty
std::shared_ptr<ExpConfig::Setup> ExpConfig::Setup::lastFound = nullptr; // default nothing found so far

template<typename T>
shared_ptr<T> ExpConfig::Get_(const THeaderInfo& header) {

    VLOG(9) << "Searching for config of type "
            << std_ext::getTypeAsString<T>();

    // check first if some manual override of the automatic
    // search is requested
    // remember if the header already contains such an information
    if(Setup::ManualName.empty() && !header.SetupName.empty()) {
        Setup::ManualName = header.SetupName;
        LOG(INFO) << "Manually set setup name " << Setup::ManualName
                  << " obtained from header " << header;
    }

    shared_ptr<Setup> config = nullptr;

    if(!Setup::ManualName.empty()) {
        config = ExpConfig::Setup::Get(Setup::ManualName);
        if(config == nullptr) {
            throw Exception(
                        std_ext::formatter()
                        << "Found no config matching name "
                        << Setup::ManualName
                        );
        }
    }
    else {
        // go to automatic search mode in all registered setups

        std::list< std::shared_ptr<Setup> > modules;
        for(auto setup_name : expconfig::SetupRegistry::GetNames()) {
            modules.emplace_back(expconfig::SetupRegistry::GetSetup(setup_name));
        }

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
        Setup::lastFound = modules.back();
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



shared_ptr<ExpConfig::Setup> ExpConfig::Setup::Get(const std::string& name)
{
    lastFound = expconfig::SetupRegistry::GetSetup(name);
    return lastFound;
}

shared_ptr<ExpConfig::Setup> ExpConfig::Setup::GetLastFound()
{
    // try if we have a name
    if(lastFound==nullptr && !ManualName.empty()) {
        Get(ManualName);
    }
    return lastFound;
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
    lastFound = nullptr;
    ManualName = "";
    expconfig::SetupRegistry::Destroy();
}

std::list<string> ExpConfig::Setup::GetNames() {
    return expconfig::SetupRegistry::GetNames();
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




