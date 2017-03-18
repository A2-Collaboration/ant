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
ExpConfig::Setup::SetupPtr ExpConfig::Setup::currentSetup = nullptr; // default nothing found so far

void ExpConfig::Setup::SetByTID(const TID& tid)
{
    if(!manualName.empty()) {
        // ignore requests to set by TID if name was specified
        return;
    }

    // go to automatic search mode in all registered setups

    std::list<SetupPtr> setups;
    for(auto setup_name : expconfig::SetupRegistry::GetNames()) {
        setups.emplace_back(expconfig::SetupRegistry::GetSetup(setup_name));
    }

    // remove the config if the config says it does not match
    setups.remove_if([&tid] (const SetupPtr& m) {
        return !m->Matches(tid);
    });

    // check if exactly one setup is left, provide helpful exception message
    if(setups.empty())
        throw Exception(std_ext::formatter() << "No setup found for TID " << tid);
    if(setups.size()>1)
        throw Exception(std_ext::formatter() << "More than one setup found for TID " << tid);

    // remember last found
    currentSetup = setups.back();

    LOG(INFO) << "Auto-detected setup with name " << currentSetup->GetName();
}

const expconfig::Setup_traits& ExpConfig::Setup::Get()
{
    if(!currentSetup)
        throw ExceptionNoSetup("No setup available. Call ExpConfig::Setup::SetBy* methods");
    return *currentSetup;
}

shared_ptr<Detector_t> ExpConfig::Setup::GetDetector(Detector_t::Type_t type)
{
    auto& config = Get();
    for(const auto& detector : config.GetDetectors()) {
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





