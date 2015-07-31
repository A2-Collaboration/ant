#include "Physics.h"
#include "base/Logger.h"

#include <iomanip>
#include <chrono>

using namespace std;
using namespace ant;

void ant::DebugPhysics::ProcessEvent(const ant::Event &event)
{
    VLOG(8) << event;
}

void ant::DebugPhysics::Finish()
{
    VLOG(8) << "Nop";
}

void ant::DebugPhysics::ShowResult()
{
    VLOG(8) << "Nop";
}


ant::Physics::Physics(const string &name):
    HistFac(name)
{}


ant::PhysicsManager::PhysicsManager() : physics()
{

}

void PhysicsManager::ReadFrom(list< unique_ptr<input::DataReader> > readers,
        long long maxevents,
        bool& running)
{

    chrono::time_point<std::chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    long long nEvents = 0;

    while(true) {
        if(!running)
            break;
        if(nEvents>=maxevents)
            break;



        nEvents++;
    }


    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    LOG(INFO) << "Processed " << nEvents << " events, speed "
              << nEvents/elapsed_seconds.count() << " Items/s";
}



void ant::PhysicsManager::ProcessEvent(const ant::Event &event)
{
    for( auto& m : physics ) {

        m->ProcessEvent(event);
    }
}

void ant::PhysicsManager::ShowResults()
{
    for(auto& p : physics) {
        p->ShowResult();
    }
}


ant::PhysicsRegistry&ant::PhysicsRegistry::get()
{
    static PhysicsRegistry instance;
    return instance;
}

std::unique_ptr<ant::Physics> ant::PhysicsRegistry::Create(const string& name)
{
    return PhysicsRegistry::get().physics_creators.at(name)();

}

void ant::PhysicsRegistry::PrintRegistry()
{
    for(auto& entry : get().physics_creators) {
        LOG(INFO) << entry.first;
    }
}


ant::PhysicsRegistration::PhysicsRegistration(ant::physics_creator c, const string& name)
{
    PhysicsRegistry::get().RegisterPhysics(c,name);
}

AUTO_REGISTER_PHYSICS(DebugPhysics, "DebugPhysics")
