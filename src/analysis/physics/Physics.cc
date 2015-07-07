#include "Physics.h"

#include "TFile.h"
#include "TDirectory.h"

#include <iostream>

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
        auto event = reader.ReadNextEvent();
        for( auto& m : physics ) {
            m->ProcessEvent(*event.get());
        }
    }
}
