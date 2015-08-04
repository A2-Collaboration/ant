#include "catch.hpp"

#include "analysis/input/goat/GoatReader.h"
#include "analysis/data/Event.h"
#include "analysis/physics/Physics.h"

#include "base/Logger.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"

#include "TH1D.h"

#include <iostream>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("InputModule", "[analysis]") {
    dotest();
}

void dotest() {

    ant::tmpfile_t tmp1;
    ant::tmpfile_t tmp2;

    WrapTFile masterFile(tmp1.filename,
                         WrapTFile::mode_t::recreate,
                         true);

    WrapTFile sndFile(tmp2.filename);
    auto h2 = sndFile.CreateInside<TH1D>("b","B",10,0,10);
    h2->Fill(3);

    auto h1 = new TH1D("a","A",10,0,10);
    h1->Fill(2);


    PhysicsManager pm;

    pm.AddPhysics<DebugPhysics>("aaaa");
}
