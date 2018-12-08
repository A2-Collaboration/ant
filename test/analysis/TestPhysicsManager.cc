#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "analysis/physics/PhysicsManager.h"
#include "analysis/input/ant/AntReader.h"
#include "analysis/input/pluto/PlutoReader.h"
#include "analysis/input/goat/GoatReader.h"

#include "analysis/utils/Uncertainties.h"
#include "analysis/utils/ParticleTools.h"
#include "analysis/utils/ParticleID.h"

#include "unpacker/Unpacker.h"
#include "reconstruct/Reconstruct.h"
#include "expconfig/ExpConfig.h"

#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"

#include "TTree.h"


#include <iostream>
#include <list>

using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest_raw();
void dotest_raw_nowrite();
void dotest_plutogeant(bool insertGoat, bool checktaggerhits = false);
void dotest_pluto(bool insertGoat);
void dotest_runall();

TEST_CASE("PhysicsManager: Raw Input", "[analysis]") {
    test::EnsureSetup();
    dotest_raw();
}

TEST_CASE("PhysicsManager: Raw Input without TEvent writing", "[analysis]") {
    test::EnsureSetup();
    dotest_raw_nowrite();
}

TEST_CASE("PhysicsManager: Pluto/Geant Input", "[analysis]") {
    test::EnsureSetup();
    dotest_plutogeant(false);
}

TEST_CASE("PhysicsManager: Pluto/Geant Input check tagger hits", "[analysis]") {
    test::EnsureSetup();
    dotest_plutogeant(false, true);
}

TEST_CASE("PhysicsManager: Pluto only Input", "[analysis]") {
    test::EnsureSetup();
    dotest_pluto(false);
}

TEST_CASE("PhysicsManager: Pluto/Geant Input with Goat", "[analysis]") {
    test::EnsureSetup();
    dotest_plutogeant(true);
}

TEST_CASE("PhysicsManager: Pluto only Input with Goat", "[analysis]") {
    test::EnsureSetup();
    dotest_pluto(true);
}

TEST_CASE("PhysicsManager: Run all physics", "[analysis]") {
    test::EnsureSetup();
    dotest_runall();
}

struct TestPhysics : Physics
{
    bool finishCalled = false;
    bool showCalled = false;
    bool nowrite    = false;
    bool checktaggerhits = false;
    unsigned seenEvents = 0;
    unsigned seenTaggerHits = 0;
    unsigned seenCandidates = 0;
    unsigned seenMCTrue = 0;
    unsigned seenTrueTargetPos = 0;
    unsigned seenReconTargetPosNaN = 0;



    TestPhysics(bool nowrite_ = false, bool checktaggerhits_ = false) :
        Physics("TestPhysics", nullptr),
        nowrite(nowrite_),
        checktaggerhits(checktaggerhits_)
    {
        HistFac.makeTH1D("test","test","test",BinSettings(10));
    }

    virtual void ProcessEvent(const TEvent& event, physics::manager_t& manager) override
    {
        seenEvents++;
        seenTaggerHits += event.Reconstructed().TaggerHits.size();
        seenCandidates += event.Reconstructed().Candidates.size();
        seenMCTrue += utils::ParticleTypeList::Make(event.MCTrue().ParticleTree).GetAll().size();
        // make sure it's non-zero and not nan only for MCTrue
        seenTrueTargetPos += event.MCTrue().Target.Vertex.z < -1;
        seenReconTargetPosNaN += std::isnan(event.Reconstructed().Target.Vertex.z);
        // request to save every third event
        if(!nowrite && seenEvents % 3 == 0)
            manager.SaveEvent();
        if(checktaggerhits) {
            auto& mctrue = event.MCTrue().TaggerHits;
            auto& recon = event.Reconstructed().TaggerHits;
            REQUIRE_FALSE(mctrue.empty());
            REQUIRE_FALSE(recon.empty());
            auto it_mctrue = mctrue.begin();
            auto it_recon = recon.begin();
            while(it_mctrue != mctrue.end() && it_recon != recon.end()) {
                REQUIRE(it_mctrue->Channel == it_recon->Channel);
                REQUIRE(it_mctrue->PhotonEnergy == Approx(it_recon->PhotonEnergy).epsilon(0.005));
                ++it_mctrue;
                ++it_recon;
            }
        }
    }
    virtual void Finish() override
    {
        finishCalled = true;
    }
    virtual void ShowResult() override
    {
        showCalled = true;
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
    const unsigned expectedEvents = 221;

    tmpfile_t tmpfile;

    // write out some file
    {
        WrapTFileOutput outfile(tmpfile.filename, true);
        PhysicsManagerTester pm;
        pm.AddPhysics<TestPhysics>();

        // make some meaningful input for the physics manager
        auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
        auto reconstruct = std_ext::make_unique<Reconstruct>();
        list< unique_ptr<analysis::input::DataReader> > readers;
        readers.emplace_back(std_ext::make_unique<input::AntReader>(nullptr, move(unpacker), move(reconstruct)));
        pm.ReadFrom(move(readers), numeric_limits<long long>::max());

        const std::uint32_t timestamp = 1408221194;
        REQUIRE(pm.GetProcessedTIDRange() == interval<TID>(TID(timestamp, 0u), TID(timestamp, expectedEvents-1)) );

        std::shared_ptr<TestPhysics> physics = pm.GetTestPhysicsModule();

        REQUIRE(physics->finishCalled);
        REQUIRE_FALSE(physics->showCalled);

        REQUIRE(physics->seenEvents == expectedEvents);
        REQUIRE(physics->seenCandidates == 864);
        REQUIRE(physics->seenTrueTargetPos == 0);
        REQUIRE(physics->seenReconTargetPosNaN == expectedEvents);

        // quick check if TTree was there...
        auto tree = outfile.GetSharedClone<TTree>("treeEvents");
        REQUIRE(tree != nullptr);
        REQUIRE(tree->GetEntries() == expectedEvents/3);
    }

    // read in file with AntReader
    {
        auto inputfiles = make_shared<WrapTFileInput>(tmpfile.filename);

        PhysicsManagerTester pm;
        pm.AddPhysics<TestPhysics>();

        // make some meaningful input for the physics manager

        auto reconstruct = std_ext::make_unique<Reconstruct>();
        list< unique_ptr<analysis::input::DataReader> > readers;

        readers.emplace_back(std_ext::make_unique<input::AntReader>(inputfiles, nullptr, move(reconstruct)));
        pm.ReadFrom(move(readers), numeric_limits<long long>::max());

        // note that we actually requested every third event to be saved in the physics class
        const std::uint32_t timestamp = 1408221194;
        REQUIRE(pm.GetProcessedTIDRange() == interval<TID>(TID(timestamp, 2u), TID(timestamp, 3*unsigned((expectedEvents-1)/3)-1)));

        std::shared_ptr<TestPhysics> physics = pm.GetTestPhysicsModule();

        REQUIRE(physics->seenEvents == expectedEvents/3);
        // make sure the reconstruction wasn't applied twice!
        REQUIRE(physics->seenCandidates == 286);

    }

    // read in file with AntReader without reconstruction
    {
        auto inputfiles = make_shared<WrapTFileInput>(tmpfile.filename);

        PhysicsManagerTester pm;
        pm.AddPhysics<TestPhysics>();

        // make some meaningful input for the physics manager

        list< unique_ptr<analysis::input::DataReader> > readers;
        readers.emplace_back(std_ext::make_unique<input::AntReader>(inputfiles, nullptr, nullptr));
        pm.ReadFrom(move(readers), numeric_limits<long long>::max());

        // note that we actually requested every third event to be saved in the physics class
        const std::uint32_t timestamp = 1408221194;
        REQUIRE(pm.GetProcessedTIDRange() == interval<TID>(TID(timestamp, 2u), TID(timestamp, 3*unsigned((expectedEvents-1)/3)-1)));

        std::shared_ptr<TestPhysics> physics = pm.GetTestPhysicsModule();

        REQUIRE(physics->seenEvents == expectedEvents/3);
        REQUIRE(physics->seenCandidates == 286);

    }

}

void dotest_raw_nowrite()
{
    tmpfile_t tmpfile;
    WrapTFileOutput outfile(tmpfile.filename, true);

    PhysicsManagerTester pm;
    pm.AddPhysics<TestPhysics>(true);

    // make some meaningful input for the physics manager
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
    auto reconstruct = std_ext::make_unique<Reconstruct>();
    list< unique_ptr<analysis::input::DataReader> > readers;
    readers.emplace_back(std_ext::make_unique<input::AntReader>(nullptr, move(unpacker), move(reconstruct)));
    pm.ReadFrom(move(readers), numeric_limits<long long>::max());

    const std::uint32_t timestamp = 1408221194;
    const unsigned expectedEvents = 221;
    REQUIRE(pm.GetProcessedTIDRange() == interval<TID>(TID(timestamp, 0u), TID(timestamp, expectedEvents-1)));

    std::shared_ptr<TestPhysics> physics = pm.GetTestPhysicsModule();

    REQUIRE(physics->finishCalled);
    REQUIRE_FALSE(physics->showCalled);

    CHECK(physics->seenEvents == expectedEvents);
    CHECK(physics->seenTaggerHits == 6272);
    CHECK(physics->seenCandidates == 864);

    // the PhysicsManager should not create a TTree...
    REQUIRE(outfile.GetSharedClone<TTree>("treeEvents") == nullptr);
}

void dotest_plutogeant(bool insertGoat, bool checktaggerhits)
{
    tmpfile_t tmpfile;
    WrapTFileOutput outfile(tmpfile.filename, true);

    PhysicsManagerTester pm;
    pm.AddPhysics<TestPhysics>(false, checktaggerhits);

    // make some meaningful input for the physics manager
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Geant_with_TID.root");
    auto reconstruct = std_ext::make_unique<Reconstruct>();
    list< unique_ptr<analysis::input::DataReader> > readers;
    readers.emplace_back(std_ext::make_unique<input::AntReader>(nullptr, move(unpacker), move(reconstruct)));

    auto rootfiles = std::make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_with_TID.root");
    readers.push_back(std_ext::make_unique<analysis::input::PlutoReader>(rootfiles));

    // should be auto-detected that Goat reader is not needed
    if(insertGoat)
        readers.push_back(std_ext::make_unique<analysis::input::GoatReader>(rootfiles));

    REQUIRE_NOTHROW(pm.ReadFrom(move(readers), numeric_limits<long long>::max()));

    std::shared_ptr<TestPhysics> physics = pm.GetTestPhysicsModule();

    CHECK(physics->seenEvents == 100);
    CHECK(physics->seenTaggerHits == 100);
    CHECK(physics->seenCandidates == 309);
    CHECK(physics->seenMCTrue == 300);
    CHECK(physics->seenTrueTargetPos == 46);
    CHECK(physics->seenReconTargetPosNaN == 100);
}

void dotest_pluto(bool insertGoat)
{
    tmpfile_t tmpfile;
    WrapTFileOutput outfile(tmpfile.filename, true);

    PhysicsManagerTester pm;
    pm.AddPhysics<TestPhysics>();

    // make some meaningful input for the physics manager
    list< unique_ptr<analysis::input::DataReader> > readers;
    auto rootfiles = std::make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_with_TID.root");
    readers.push_back(std_ext::make_unique<analysis::input::PlutoReader>(rootfiles));

    // should be auto-detected that Goat reader is not needed
    if(insertGoat)
        readers.push_back(std_ext::make_unique<analysis::input::GoatReader>(rootfiles));

    REQUIRE_NOTHROW(pm.ReadFrom(move(readers), numeric_limits<long long>::max()));

    std::shared_ptr<TestPhysics> physics = pm.GetTestPhysicsModule();

    CHECK(physics->seenEvents == 100);
    CHECK(physics->seenTaggerHits == 0); // recon only
    CHECK(physics->seenCandidates == 0); // recon only
    CHECK(physics->seenMCTrue == 300);
    CHECK(physics->seenTrueTargetPos == 0); // pluto does not know anything about target
    CHECK(physics->seenReconTargetPosNaN == 100);
}

void dotest_runall() {

    // ensure simple particle ID, some physics classes need it
    utils::ParticleID::SetDefault(std_ext::make_unique<utils::SimpleParticleID>());

    for(auto name : PhysicsRegistry::GetList()) {

        // Exclude known physics classes which implement method calls which will segfault when using ROOT6
        /// \todo Check if the bug reported <a href="https://sft.its.cern.ch/jira/browse/ROOT-9450">here</a> is fixed and remove exclusion
        /// \bug Certain physics classes will probably crash due to <a href="https://sft.its.cern.ch/jira/browse/ROOT-9450">ROOT-9450</a>
        #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
        bool badPhysics = false;
        for (const string& physics : {"EtapDalitz", "EtapDalitzMC", "EtapOmegaG", "Omega_EpEm", "sigmaPlus", "singlePi0", "triplePi0"})
            if (!physics.compare(name))
                badPhysics = true;
        if (badPhysics)
            continue;
        #endif

        PhysicsManager pm;
        INFO(name);
        try {
            pm.AddPhysics(PhysicsRegistry::Create(name));
        }
        catch(ExpConfig::ExceptionNoDetector&) {
            // ignore silently if test setup did not provide detector
            continue;
        }
        catch(utils::UncertaintyModel::Exception&) {
            // ignore silently if class cannot load uncertainty model
            continue;
        }
        catch(Physics::ExceptionOptionNeeded&) {
            // ignore silently if class needs user option
            continue;
        }
        catch(const std::exception& e) { // use reference to prevent object slicing!
            FAIL(string("Unexpected exception while creating physics class: ")+e.what());
        }
        catch(...) {
            FAIL("Something weird was thrown.");
        }
        auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_twoscalerblocks.dat.xz");
        auto reconstruct = std_ext::make_unique<Reconstruct>();
        list< unique_ptr<analysis::input::DataReader> > readers;
        readers.emplace_back(std_ext::make_unique<input::AntReader>(nullptr, move(unpacker), move(reconstruct)));

        REQUIRE_NOTHROW(pm.ReadFrom(move(readers), 20));
    }

}

