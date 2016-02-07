#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "analysis/physics/PhysicsManager.h"
#include "analysis/input/ant/AntReader.h"
#include "analysis/input/pluto/PlutoReader.h"

#include "analysis/slowcontrol/SlowControlVariables.h"

#include "unpacker/Unpacker.h"
#include "reconstruct/Reconstruct.h"

#include "base/Logger.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"

#include "TTree.h"


#include <iostream>
#include <list>

using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest();

TEST_CASE("SlowControlManager", "[analysis]") {
    test::EnsureSetup();
    dotest();
}


struct TestPhysics : Physics
{


    TestPhysics() :
        Physics("TestPhysics", nullptr)
    {
        slowcontrol::Variables::TaggerScalers->Request();
    }

    virtual void ProcessEvent(const TEvent& event, physics::manager_t& manager) override
    {
        REQUIRE_FALSE(event.SavedForSlowControls);
        auto taggerscalers = slowcontrol::Variables::TaggerScalers->Get();
        REQUIRE(taggerscalers.size() == 47);
        if(event.Reconstructed->SlowControls.empty())
            manager.SaveEvent();
    }
};

void dotest()
{

    tmpfile_t stage1;
    tmpfile_t stage2;

    // write out some file
    {
        WrapTFileOutput outfile(stage1.filename, WrapTFileOutput::mode_t::recreate, true);
        PhysicsManager pm;
        pm.AddPhysics<TestPhysics>();

        // make some meaningful input for the physics manager
        auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_twoscalerblocks.dat.xz");
        auto reconstruct = std_ext::make_unique<Reconstruct>();
        list< unique_ptr<analysis::input::DataReader> > readers;
        readers.emplace_back(std_ext::make_unique<input::AntReader>(nullptr, move(unpacker), move(reconstruct)));
        pm.ReadFrom(move(readers), numeric_limits<long long>::max());

        // quick check if TTree was there...
        auto tree = outfile.GetSharedClone<TTree>("treeEvents");
        REQUIRE(tree != nullptr);
        REQUIRE(tree->GetEntries() == 422);
    }

    // read in file with AntReader, write to stage2
    {
        WrapTFileOutput outfile(stage2.filename, WrapTFileOutput::mode_t::recreate, true);

        auto inputfiles = make_shared<WrapTFileInput>(stage1.filename);

        PhysicsManager pm;
        pm.AddPhysics<TestPhysics>();

        // make some meaningful input for the physics manager

        auto reconstruct = std_ext::make_unique<Reconstruct>();
        list< unique_ptr<analysis::input::DataReader> > readers;

        readers.emplace_back(std_ext::make_unique<input::AntReader>(inputfiles, nullptr, move(reconstruct)));
        pm.ReadFrom(move(readers), numeric_limits<long long>::max());

        // quick check if TTree was there...
        auto tree = outfile.GetSharedClone<TTree>("treeEvents");
        REQUIRE(tree != nullptr);
        REQUIRE(tree->GetEntries() == 422);
    }

    // read in file with AntReader
    {
        auto inputfiles = make_shared<WrapTFileInput>(stage2.filename);

        PhysicsManager pm;
        pm.AddPhysics<TestPhysics>();

        // make some meaningful input for the physics manager

        list< unique_ptr<analysis::input::DataReader> > readers;
        readers.emplace_back(std_ext::make_unique<input::AntReader>(inputfiles, nullptr, nullptr));
        pm.ReadFrom(move(readers), numeric_limits<long long>::max());
    }

}
