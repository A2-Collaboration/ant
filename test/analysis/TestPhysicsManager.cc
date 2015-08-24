#include "catch.hpp"

#include "analysis/input/goat/GoatReader.h"
#include "analysis/data/Event.h"
#include "analysis/physics/PhysicsManager.h"
#include "analysis/physics/common/DebugPhysics.h"

#include "base/Logger.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"

#include "TH1D.h"

#include <iostream>
#include <list>

using namespace std;
using namespace ant;

void dotest_write(const string& filename1,const string& filename2);
void dotest_read(const string& filename1,const string& filename2);

TEST_CASE("InputModule", "[analysis]") {
    tmpfile_t tmp1;
    tmpfile_t tmp2;
    dotest_write(tmp1.filename,tmp2.filename);

    dotest_read(tmp1.filename,tmp2.filename);
}

void dotest_write(const string& filename1, const string& filename2)
{
    WrapTFileOutput masterFile(filename1,
                         WrapTFileOutput::mode_t::recreate,
                         true);

    WrapTFileOutput sndFile(filename2);
    auto h2 = sndFile.CreateInside<TH1D>("b","B",10,0,10);
    auto h3 = sndFile.CreateInside<TH1D>("c","C",10,0,10);
    h2->Fill(3);

    auto h1 = new TH1D("a","A",10,0,10);
    h1->Fill(2);

    analysis::PhysicsManager pm;

    pm.AddPhysics<analysis::physics::DebugPhysics>();
}

void dotest_read(const string& filename1, const string& filename2)
{
    WrapTFileInput masterRead(filename1);
    WrapTFileInput sndRead(filename2);

    REQUIRE(masterRead.GetListOf<TH1D>().size() == 1);
    REQUIRE(sndRead.GetListOf<TH1D>().size() == 2);
}
