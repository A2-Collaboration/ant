#include "catch.hpp"

#include "analysis/input/goat/GoatReader.h"
#include "analysis/data/Event.h"
#include "analysis/physics/Physics.h"

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
    WrapTFile masterFile(filename1,
                         WrapTFile::mode_t::recreate,
                         true);

    WrapTFile sndFile(filename2);
    auto h2 = sndFile.CreateInside<TH1D>("b","B",10,0,10);
    auto h3 = sndFile.CreateInside<TH1D>("c","C",10,0,10);
    h2->Fill(3);

    auto h1 = new TH1D("a","A",10,0,10);
    h1->Fill(2);

    PhysicsManager pm;

    pm.AddPhysics<DebugPhysics>("aaaa");
}

void dotest_read(const string& filename1, const string& filename2)
{
    WrapTFile masterRead(filename1,
                         WrapTFile::mode_t::read);
    WrapTFile sndRead(filename2,
                      WrapTFile::mode_t::read);

    REQUIRE(masterRead.GetListOf<TH1D>().size() == 1);
    REQUIRE(sndRead.GetListOf<TH1D>().size() == 2);
}
