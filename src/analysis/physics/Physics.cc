#include "Physics.h"

#include "TFile.h"
#include "TDirectory.h"
#include "base/Logger.h"

#include <iomanip>

using namespace std;

void ant::DebugPhysics::ProcessEvent(const ant::Event &event)
{
    cout << event << endl;
}

void ant::DebugPhysics::Finish()
{
    cout << "DebugPhysics Done." << endl;
}

void ant::DebugPhysics::ShowResult()
{
    cout << "DebugPhysics ShowResults." << endl;
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
        for( auto& m : physics ) {

            m->ProcessEvent(*event.get());
        }

        const auto i = reader.EventsRead();
        if( i % 10000 == 0) {
            const auto nevents = reader.TotalEvents();
            if(nevents>0)
                VLOG(3) << "Events processed: " << i << " (" << std::fixed << std::setprecision(2) << (double(i)/double(nevents)*100.0) << "%)";
            else
                VLOG(3) << "Events processed: " << i;
        }
    }
}
