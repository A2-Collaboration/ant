#include "catch.hpp"
#include "catch_config.h"

#include "analysis/data/Event.h"
#include "analysis/physics/PhysicsManager.h"
#include "analysis/input/ant/AntReader.h"
#include "analysis/input/pluto/PlutoReader.h"

#include "unpacker/Unpacker.h"
#include "reconstruct/Reconstruct.h"
#include "tree/TAntHeader.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"



#include <iostream>
#include <list>

using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest_raw();
void dotest_plutogeant();


TEST_CASE("PhysicsManager: Raw Input", "[analysis]") {
    dotest_raw();
}

TEST_CASE("PhysicsManager: Pluto/Geant Input", "[analysis]") {
    dotest_plutogeant();
}

struct TestPhysics : Physics
{
    bool finishCalled = false;
    bool initCalled = false;
    bool showCalled = false;
    unsigned seenEvents = 0;
    unsigned seenCandidates = 0;
    unsigned seenMCTrue = 0;


    TestPhysics() : Physics("TestPhysics") {}

    virtual void ProcessEvent(const data::Event& event) override
    {
        seenEvents++;
        seenCandidates += event.Reconstructed().Candidates().size();
        seenMCTrue += event.MCTrue().Particles().GetAll().size();
    }
    virtual void Finish() override
    {
        finishCalled = true;
    }
    virtual void ShowResult() override
    {
        showCalled = true;
    }
    virtual void Initialize(data::Slowcontrol&) override
    {
        initCalled = true;
    }
};

struct PhysicsManagerTester : PhysicsManager
{
    using PhysicsManager::PhysicsManager;

    shared_ptr<TestPhysics> GetTestPhysicsModule() {
        // bit ugly to obtain the physics module back
        auto module = move(physics.back());
        physics.pop_back();
        return dynamic_pointer_cast<TestPhysics, Physics>(
                    std::shared_ptr<Physics>(module.release()));
    }
};

void dotest_raw()
{
    PhysicsManagerTester pm;
    pm.AddPhysics<TestPhysics>();

    // make some meaningful input for the physics manager
    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
    auto reconstruct = std_ext::make_unique<Reconstruct>();
    list< unique_ptr<analysis::input::DataReader> > readers;
    readers.emplace_back(std_ext::make_unique<input::AntReader>(move(unpacker), move(reconstruct)));
    bool running = true;
    TAntHeader header;
    pm.ReadFrom(move(readers), numeric_limits<long long>::max(), running, header);

    const std::uint32_t timestamp = 1408221194;
    const unsigned expectedEvents = 221;
    REQUIRE(header.FirstID == TID(timestamp, 0u));
    REQUIRE(header.LastID == TID(timestamp, expectedEvents-1));

    std::shared_ptr<TestPhysics> physics = pm.GetTestPhysicsModule();

    REQUIRE(physics->finishCalled);
    REQUIRE_FALSE(physics->showCalled);
    REQUIRE(physics->initCalled);

    REQUIRE(physics->seenEvents == expectedEvents);
    REQUIRE(physics->seenCandidates == 822);
}

void dotest_plutogeant()
{
    PhysicsManagerTester pm;
    pm.AddPhysics<TestPhysics>();

    // make some meaningful input for the physics manager
    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Geant_with_TID.root");
    auto reconstruct = std_ext::make_unique<Reconstruct>();
    list< unique_ptr<analysis::input::DataReader> > readers;
    readers.emplace_back(std_ext::make_unique<input::AntReader>(move(unpacker), move(reconstruct)));

    auto plutofile = std::make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_with_TID.root");
    readers.push_back(std_ext::make_unique<analysis::input::PlutoReader>(plutofile, std::shared_ptr<TaggerDetector_t>(nullptr)));

    bool running = true;
    TAntHeader header;
    pm.ReadFrom(move(readers), numeric_limits<long long>::max(), running, header);

    std::shared_ptr<TestPhysics> physics = pm.GetTestPhysicsModule();

    REQUIRE(physics->seenEvents == 10);
    REQUIRE(physics->seenCandidates == 10);
    REQUIRE(physics->seenMCTrue == 10);

}
