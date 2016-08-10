#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "analysis/slowcontrol/SlowControlManager.h"

#include "analysis/physics/PhysicsManager.h"
#include "analysis/input/ant/AntReader.h"

#include "analysis/slowcontrol/SlowControlProcessors.h"
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

void dotest_ScalerBlobs();
void dotest_FakeReader();

//TEST_CASE("SlowControlManager: Two scaler blob", "[analysis]") {
//    test::EnsureSetup();
//    dotest_ScalerBlobs();
//}

TEST_CASE("SlowControlManager: FakeReader", "[analysis]") {
    dotest_FakeReader();
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
        auto taggerscalers = slowcontrol::Variables::TaggerScalers->Get();
        REQUIRE(taggerscalers.size() == 47);
        if(event.Reconstructed().SlowControls.empty())
            manager.SaveEvent();
    }
};

void dotest_ScalerBlobs()
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
        REQUIRE(tree->GetEntries() == 212);
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
        REQUIRE(tree->GetEntries() == 212);
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

// define some test processors

unsigned maxEvents = 16; // TIDs from 0x0 to 0xf, good for debugging

struct TestProcessor : slowcontrol::Processor {
    unsigned nProcessed = 0;
    unsigned nPopped = 0;
    unsigned nCompleted = 0;
    virtual return_t ProcessEventData(const TEventData&, physics::manager_t&) override {
        nProcessed++;
        return {};
    }
    virtual void PopQueue() override {
        nPopped++;
    }
    virtual ~TestProcessor() {
        REQUIRE(nProcessed==maxEvents);
        REQUIRE(nCompleted>0);
    }
};

struct TestProcessor1 : TestProcessor {

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override
    {
        TestProcessor::ProcessEventData(recon, manager);
        if(recon.ID.Timestamp < 0x2)
            return return_t::Skip;
        if(recon.ID.Timestamp == 0x6 || recon.ID.Timestamp == 0x9) {
            manager.SaveEvent();
            nCompleted++;
            return return_t::Complete;
        }
        return return_t::Buffer;
    }
    virtual ~TestProcessor1() {
        REQUIRE(nPopped==0);
    }
};

struct TestProcessor2 : TestProcessor {

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override
    {
        TestProcessor::ProcessEventData(recon, manager);
        if(recon.ID.Timestamp < 0x5)
            return return_t::Skip;
        if(recon.ID.Timestamp == 0x6 || recon.ID.Timestamp == 0xa) {
            manager.SaveEvent();
            nCompleted++;
            return return_t::Complete;
        }
        return return_t::Buffer;
    }
};

struct TestProcessor3 : TestProcessor {

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override
    {
        TestProcessor::ProcessEventData(recon, manager);
        if(recon.ID.Timestamp == 0x0) {
            manager.SaveEvent();
            nCompleted++;
            return return_t::Complete;
        }
        return return_t::Process;
    }
};

struct TestProcessor4 : TestProcessor {

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override
    {
        TestProcessor::ProcessEventData(recon, manager);
        if(recon.ID.Timestamp < 0x2) {
            manager.SaveEvent();
            return return_t::Skip;
        }
        else if(recon.ID.Timestamp == 0x2) {
            nCompleted++;
            return return_t::Complete;
        }
        return return_t::Process;
    }
};

struct TestSlowControlManager : SlowControlManager {
    TestSlowControlManager() : SlowControlManager() {
        // previous tests might have requested static slowcontrol variables
        // and the default ctor searches for it...
        slowcontrol.clear();
        // we add our own test processors
        AddProcessor(make_shared<TestProcessor1>());
        AddProcessor(make_shared<TestProcessor2>());
        AddProcessor(make_shared<TestProcessor3>());
        AddProcessor(make_shared<TestProcessor4>());
        REQUIRE(slowcontrol.size() == 4);
    }
};

void dotest_FakeReader() {

    TestSlowControlManager scm;

    // this is basically how PhysicsManager drives the SlowControlManager
    unsigned nEventsRead=0;
    unsigned nEventsPopped = 0;
    while(nEventsRead<maxEvents) {
        while(nEventsRead<maxEvents) {
            TEvent event;
            event.MakeReconstructed(nEventsRead);
            ++nEventsRead; // break might occur, but increase i
            if(scm.ProcessEvent(move(event)))
                break;
        }

        while(auto event = scm.PopEvent()) {
            REQUIRE(event.Event.HasReconstructed());
            nEventsPopped++;
        }
    }

    CHECK(nEventsPopped == maxEvents);
    CHECK(nEventsRead == maxEvents);

}
