#include "AntPhysics.h"

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
