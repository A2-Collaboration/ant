#include "catch.hpp"

#include "analysis/input/goat/GoatReader.h"
#include "analysis/data/Event.h"
#include "analysis/OutputManager.h"
#include "analysis/physics/Physics.h"

#include "base/Logger.h"
#include "base/tmpfile_t.h"

#include "TH1D.h"

#include <iostream>

using namespace std;
using namespace ant::output;
using namespace ant;

void dotest();

TEST_CASE("InputModule", "[analysis]") {
    dotest();
}

void dotest() {

    OutputManager om;

    ant::tmpfile_t tmp1;
    ant::tmpfile_t tmp2;

    om.SetNewOutput(tmp1.filename);
    auto h1 = new TH1D("a","A",10,0,10);
    h1->Fill(2);

    om.SetNewOutput(tmp2.filename);
    auto h2 = new TH1D("b","B",10,0,10);
    h2->Fill(3);

    PhysicsManager pm;

    pm.AddPhysics<DebugPhysics>("aaaa");
}
