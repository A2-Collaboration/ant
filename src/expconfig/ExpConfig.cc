#include "ExpConfig.h"

#include "setups/Setup.h"

#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"

#include "tree/TID.h"
#include "base/Logger.h"

#include <type_traits>
#include <list>
#include <iostream>

using namespace std;
using namespace ant;

std::string ExpConfig::Setup::manualName = ""; // default empty
ExpConfig::SetupPtr ExpConfig::Setup::currentSetup = nullptr; // default nothing found so far

void ExpConfig::Setup::SetByTID(const TID& tid)
{
    if(!manualName.empty()) {
        // ignore requests to set by TID if name was specified
        return;
    }

    // go to automatic search mode in all registered setups

    std::list<ExpConfig::SetupPtr> modules;
    for(auto setup_name : expconfig::SetupRegistry::GetNames()) {
        modules.emplace_back(expconfig::SetupRegistry::GetSetup(setup_name));
    }

    // remove the config if the config says it does not match
    modules.remove_if([&tid] (const ExpConfig::SetupPtr& m) {
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

    // remember last found
    currentSetup = modules.back();

    LOG(INFO) << "Auto-detected setup with name " << currentSetup->GetName();
}

ExpConfig::SetupPtr ExpConfig::Setup::Get()
{
    if(!currentSetup)
        throw ExceptionNoSetup("No setup available. Call ExpConfig::Setup::SetBy* methods");
    return currentSetup;
}

shared_ptr<Detector_t> ExpConfig::Setup::GetDetector(Detector_t::Type_t type)
{
    auto config = Get();
    if(config == nullptr)
        throw ExceptionNoSetup("Could not find setup to search for required detector");
    for(const auto& detector : config->GetDetectors()) {
        if(detector->Type == type)
            return detector;
    }
    throw ExceptionNoDetector("Could not find detector in given setup");
}

void ExpConfig::Setup::SetByName(const string& setupname)
{
    manualName = setupname;
    currentSetup = expconfig::SetupRegistry::GetSetup(setupname);
    if(!currentSetup)
        throw ExceptionNoSetup("No setup found in registry with name "+setupname);
}

void ExpConfig::Setup::Cleanup()
{
    currentSetup = nullptr;
    manualName = "";
    expconfig::SetupRegistry::Cleanup();
}

std::list<string> ExpConfig::Setup::GetNames() {
    return expconfig::SetupRegistry::GetNames();
}





