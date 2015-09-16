#include "catch.hpp"

#include "analysis/data/Event.h"
#include "analysis/physics/PhysicsManager.h"

#include "base/Logger.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"


#include <iostream>
#include <list>

using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest();

TEST_CASE("PhysicsManager", "[analysis]") {
    dotest();
}

struct TestPhysics : Physics {
    bool finishCalled = false;


    TestPhysics() : Physics("TestPhysics") {}

    virtual void ProcessEvent(const data::Event& event) override
    {

    }
    virtual void Finish() override
    {
        finishCalled = true;
    }
    virtual void ShowResult() override
    {

    }
    virtual void Initialize(data::Slowcontrol&) override
    {

    }
};

void dotest()
{
    tmpfile_t tmp;
    auto outfile = std_ext::make_unique<WrapTFileOutput>(tmp.filename,
                         WrapTFileOutput::mode_t::recreate,
                         true);

    analysis::PhysicsManager pm;
    auto physics = std_ext::make_unique<TestPhysics>();
    pm.AddPhysics(move(physics));




}
