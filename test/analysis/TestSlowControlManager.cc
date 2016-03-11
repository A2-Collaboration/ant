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

void dotest_ScalerBlobs();
void dotest_FakeReader();

TEST_CASE("SlowControlManager: Two scaler blob", "[analysis]") {
    test::EnsureSetup();
    dotest_ScalerBlobs();
}

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
        REQUIRE_FALSE(event.SavedForSlowControls);
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

class SLCFakeReader : public ant::analysis::input::DataReader {
public:
    using pattern_t = std::list<std::list<TSlowControl>>;
protected:
    const pattern_t pattern;
    pattern_t::const_iterator i;

public:
    SLCFakeReader(const pattern_t& p):
        pattern(p),
        i(pattern.cbegin())
    {}

    bool IsSource() override { return true; }
    bool ReadNextEvent(TEvent &event) override;
    double PercentDone() const override { return 0.0; }


};

bool SLCFakeReader::ReadNextEvent(TEvent &event)
{
    if(i==pattern.cend())
        return false;

    for(const auto& s : *i) {
        event.Reconstructed().SlowControls.push_back(s);
    }

    ++i;

    return true;
}

TSlowControl makeAcqu(const std::string& name, const int value) {
    TSlowControl s(TSlowControl::Type_t::AcquScaler,TSlowControl::Validity_t::Backward,0, name, "Auto Gen");
    s.Payload_Int.push_back(TKeyValue<int64_t>(0,value));
    return s;
}

void dotest_FakeReader() {
    SLCFakeReader reader({
                             {makeAcqu("A", 0)},
                             {},
                             {},
                             {makeAcqu("A", 1)},
                             {},
                             {makeAcqu("A", 2),makeAcqu("B",1)},
                             {}
                         });

    TID tid(time(nullptr), 0, {TID::Flags_t::AdHoc});
    do {

        TEvent e(tid);
        ++tid;

        if(reader.ReadNextEvent(e)) {
            cout << e << endl;

        } else {
            break;
        }
    } while(true);

}
