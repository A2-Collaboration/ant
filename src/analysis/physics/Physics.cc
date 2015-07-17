#include "Physics.h"

#include "TFile.h"
#include "TDirectory.h"
#include "base/Logger.h"

#include <iomanip>

using namespace std;

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


ant::PhysicsManager::PhysicsManager()
{

}

void ant::PhysicsManager::ReadFrom(ant::input::DataReader &reader)
{
    while(reader.hasData()) {
        const auto event = reader.ReadNextEvent();
        ProcessEvent(event);

        const auto i = reader.EventsRead();
        if( i % 10000 == 0) {
            const auto nevents = reader.TotalEvents();
            if(nevents>0) {
                VLOG(3) << "Events processed: " << i << " (" << std::fixed << std::setprecision(2) << (double(i)/double(nevents)*100.0) << "%)";
            } else {
                VLOG(3) << "Events processed: " << i;
            }
        }
    }
    VLOG(3) << "No more data to read";
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
