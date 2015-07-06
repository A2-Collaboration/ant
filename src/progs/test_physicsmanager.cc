#include "input/GoatReader.h"
#include "Event.h"
#include "TH1D.h"

#include <iostream>
#include "OutputManager.h"
#include "Physics.h"

using namespace std;
using namespace ant::output;
using namespace ant;

int main(int argc, char** argv) {

    OutputManager om;

    om.SetNewOutput("out.root");
    TH1* h1 = new TH1D("a","A",10,0,10);

    om.SetNewOutput("out2.root");
    TH1* h2 = new TH1D("b","B",10,0,10);

    PhysicsManager pm;

    pm.AddPhysics<DebugPhysics>("aaaa");

    return 0;
}
