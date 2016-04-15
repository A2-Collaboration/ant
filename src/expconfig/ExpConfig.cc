#include "ExpConfig.h"

#include "setups/Setup.h"
#include "tree/TID.h"
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
shared_ptr<T> ExpConfig::Get_(const TID& tid) {

    VLOG(9) << "Searching for config of type "
            << std_ext::getTypeAsString<T>();

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
        modules.remove_if([&tid] (const shared_ptr<Setup>& m) {
            return !m->Matches(tid);
        });

        // check if something reasonable is left
        if(modules.empty()) {
            throw ExpConfig::Exception(std_ext::formatter()
                                       << "No setup found for TID "
                                       << tid);
        }
        if(modules.size()>1) {
            throw ExpConfig::Exception(std_ext::formatter()
                                       << "More than one setup found for TID "
                                       << tid);
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


shared_ptr<ExpConfig::Setup> ExpConfig::Setup::Get(const TID& tid)
{
     return Get_<ExpConfig::Setup>(tid);
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

shared_ptr<Detector_t> ExpConfig::Setup::GetDetector(Detector_t::Type_t type)
{
    auto config = GetLastFound();
    if(config == nullptr)
        throw ExceptionNoConfig("Could not find setup to search for required detector");
    for(const auto& detector : config->GetDetectors()) {
        if(detector->Type == type)
            return detector;
    }
    throw Exception("Could not find detector in given setup");
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


} // namespace ant




