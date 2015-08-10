#include "Physics.h"
#include "base/Logger.h"

#include "tree/TSlowControl.h"
#include "tree/TAntHeader.h"
#include "utils/ParticleID.h"

#include "TDirectory.h"

#include <iomanip>
#include <chrono>


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

void DebugPhysics::ProcessEvent(const data::Event& event)
{
    VLOG(8) << event;
}

void DebugPhysics::Finish()
{
    VLOG(8) << "Nop";
}

void DebugPhysics::ShowResult()
{
    VLOG(8) << "Nop";
}


Physics::Physics(const string &name):
    name_(name),
    HistFac(name)
{}


PhysicsManager::PhysicsManager() : physics()
{

}

void PhysicsManager::ReadFrom(
        list< unique_ptr<input::DataReader> > readers,
        long long maxevents,
        bool& running,
        TAntHeader* header)
{
    if(readers.empty())
        return;

    if(physics.empty()) {
        LOG(WARNING) << "No Analysis Instances activated. Will not analyse anything.";
    }


    // figure out what set of readers we have
    // we expect maximum one source and several amends

    unique_ptr<input::DataReader> source = nullptr;
    auto it_reader = readers.begin();
    while(it_reader != readers.end()) {
        if((*it_reader)->IsSource()) {
            if(source != nullptr) {
                LOG(ERROR) << "Found more than one source for events, stop.";
                return;
            }
            source = move(*it_reader);
            it_reader = readers.erase(it_reader);
        }
        ++it_reader;
    }



    chrono::time_point<std::chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    long long nEvents = 0;

    TID& firstEventID = header->FirstID;
    TID& lastEventID  = header->LastID;
    while(true) {
        if(!running)
            break;
        if(nEvents>=maxevents)
            break;

        data::Event event;
        TSlowControl slowcontrol;

        if(source) {
            if(!source->ReadNextEvent(event, slowcontrol)) {
                break;
            }
        }

        it_reader = readers.begin();
        while(it_reader != readers.end()) {
            if(!(*it_reader)->ReadNextEvent(event, slowcontrol))
                it_reader = readers.erase(it_reader);
            ++it_reader;
        }

        if(!source && readers.empty())
            break;

        if(nEvents==0)
            firstEventID = event.Reconstructed().TriggerInfos().EventID();

        lastEventID = event.Reconstructed().TriggerInfos().EventID();

        /// \todo make use of slowcontrol
        ProcessEvent(event);

        nEvents++;
    }

    VLOG(5) << "First EventId processed: " << firstEventID;
    VLOG(5) << "Last  EventId processed: " << lastEventID;

    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    LOG(INFO) << "Processed " << nEvents << " events, speed "
              << nEvents/elapsed_seconds.count() << " Items/s";
}



void PhysicsManager::ProcessEvent(data::Event &event)
{
    if(particleID) {

        // run particle ID for Reconstructed candidates
        for(auto cand : event.Reconstructed().Candidates()) {

            event.Reconstructed().Particles().AddParticle(
                        particleID->Process(cand)
                        );
        }
    }


    for( auto& m : physics ) {

        m->ProcessEvent(event);
    }
}

void PhysicsManager::ShowResults()
{
    for(auto& p : physics) {
        p->ShowResult();
    }
}


PhysicsRegistry& PhysicsRegistry::get()
{
    static PhysicsRegistry instance;
    return instance;
}

std::unique_ptr<Physics> PhysicsRegistry::Create(const string& name)
{
    return PhysicsRegistry::get().physics_creators.at(name)();

}

void PhysicsRegistry::PrintRegistry()
{
    for(auto& entry : get().physics_creators) {
        LOG(INFO) << entry.first;
    }
}


PhysicsRegistration::PhysicsRegistration(physics_creator c, const string& name)
{
    PhysicsRegistry::get().RegisterPhysics(c,name);
}

AUTO_REGISTER_PHYSICS(DebugPhysics, "DebugPhysics")
