#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"
#include "stealer.h"

#include "analysis/slowcontrol/SlowControlManager.h"

#include "analysis/physics/PhysicsManager.h"
#include "analysis/input/ant/AntReader.h"

#include "analysis/slowcontrol/SlowControlProcessors.h"
#include "analysis/slowcontrol/SlowControlVariables.h"

#include "unpacker/Unpacker.h"
#include "reconstruct/Reconstruct.h"

#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"
#include "base/std_ext/container.h"

#include "TTree.h"


#include <iostream>
#include <list>
#include <queue>

using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest_ScalerBlobs();
void dotest_FakeReader();

TEST_CASE("SlowControlManager: Two scaler blob", "[analysis]") {
    test::EnsureSetup();
    dotest_ScalerBlobs();
}

struct result_t {
    unsigned nEventsRead = 0;
    unsigned nEventsPopped = 0;
    unsigned nContextSwitched = 0;
    unsigned nEventsSkipped = 0;
    unsigned nEventsSavedForSC = 0;
};

result_t run_TestSlowControlManager(const vector<unsigned>& enabled);

TEST_CASE("SlowControlManager: Processors {1}", "[analysis]") {
    auto r = run_TestSlowControlManager({1});
    CHECK(r.nEventsPopped == 15);
    CHECK(r.nContextSwitched == 3);
    CHECK(r.nEventsSkipped == 1);
    CHECK(r.nEventsSavedForSC == 4);
}

TEST_CASE("SlowControlManager: Processors {2}", "[analysis]") {
    auto r = run_TestSlowControlManager({2});
    CHECK(r.nEventsPopped == 14);
    CHECK(r.nContextSwitched == 3);
    CHECK(r.nEventsSkipped == 1);
    CHECK(r.nEventsSavedForSC == 4);
}

TEST_CASE("SlowControlManager: Processors {3}", "[analysis]") {
    auto r = run_TestSlowControlManager({3});
    CHECK(r.nEventsPopped == 16);
    CHECK(r.nContextSwitched == 15);
    CHECK(r.nEventsSkipped == 0);
    CHECK(r.nEventsSavedForSC == 3);
}

TEST_CASE("SlowControlManager: Processors {4}", "[analysis]") {
    auto r = run_TestSlowControlManager({4});
    CHECK(r.nEventsPopped == 16);
    CHECK(r.nContextSwitched == 13);
    CHECK(r.nEventsSkipped == 2);
    CHECK(r.nEventsSavedForSC == 2);
}

TEST_CASE("SlowControlManager: Processors {1,2}", "[analysis]") {
    auto r = run_TestSlowControlManager({1,2});
    CHECK(r.nEventsPopped == 15);
    CHECK(r.nContextSwitched == 3);
    CHECK(r.nEventsSkipped == 2);
    CHECK(r.nEventsSavedForSC == 6);
}

TEST_CASE("SlowControlManager: Processors {3,4}", "[analysis]") {
    auto r = run_TestSlowControlManager({3,4});
    CHECK(r.nEventsPopped == 16);
    CHECK(r.nContextSwitched == 13);
    CHECK(r.nEventsSkipped == 2);
    CHECK(r.nEventsSavedForSC == 4);
}

TEST_CASE("SlowControlManager: Processors {1,4}", "[analysis]") {
    auto r = run_TestSlowControlManager({1,4});
    CHECK(r.nEventsPopped == 16);
    CHECK(r.nContextSwitched == 3);
    CHECK(r.nEventsSkipped == 2);
    CHECK(r.nEventsSavedForSC == 5);
}

TEST_CASE("SlowControlManager: Processors {2,3}", "[analysis]") {
    auto r = run_TestSlowControlManager({2,3});
    CHECK(r.nEventsPopped == 15);
    CHECK(r.nContextSwitched == 3);
    CHECK(r.nEventsSkipped == 2);
    CHECK(r.nEventsSavedForSC == 6);
}

TEST_CASE("SlowControlManager: Processors {1,2,3,4}", "[analysis]") {
    auto r = run_TestSlowControlManager({1,2,3,4});
    CHECK(r.nEventsPopped == 16);
    CHECK(r.nContextSwitched == 3);
    CHECK(r.nEventsSkipped == 3);
    CHECK(r.nEventsSavedForSC == 8);
}

// see https://github.com/zjx20/stealer for STEALER usage

STEALER(stealer_Variable_t, slowcontrol::Variable,
        STEAL_CONST_METHOD(list<std::shared_ptr<slowcontrol::Processor>>, GetNeededProcessors)
);

using AcquProcessor_t = slowcontrol::processor::AcquScalerVector;
STEALER(stealer_AcquProcessor_t, AcquProcessor_t,
        STEAL_FIELD(bool, firstScalerSeen),
        STEAL_FIELD(std::queue<AcquProcessor_t::value_t>, queue),
);

void reset_acquprocessors(std::shared_ptr<const slowcontrol::Variable> var) {
    stealer_Variable_t var_(*const_pointer_cast<slowcontrol::Variable>(var));
    for(auto& proc : var_.GetNeededProcessors()) {
        auto acquproc = dynamic_pointer_cast<AcquProcessor_t, slowcontrol::Processor>(proc);
        if(!acquproc)
            continue;
        stealer_AcquProcessor_t proc_(*acquproc);
        proc_.firstScalerSeen = false;
        proc_.queue = std::queue<AcquProcessor_t::value_t>();
    }
}

struct TestPhysics : Physics
{
    unsigned nChanged = 0;
    bool firstEvent = true;

    TestPhysics() :
        Physics("TestPhysics", nullptr)
    {
        // resetting important since we run SlowControlManager several times
        reset_acquprocessors(slowcontrol::Variables::Trigger);
        slowcontrol::Variables::Trigger->Request();
        reset_acquprocessors(slowcontrol::Variables::TaggerScalers);
        slowcontrol::Variables::TaggerScalers->Request();

    }

    virtual void ProcessEvent(const TEvent& event, physics::manager_t& manager) override
    {
        if(slowcontrol::Variables::TaggerScalers->HasChanged()) {
            CHECK(firstEvent);
            nChanged++;
            CHECK(event.Reconstructed().ID.Lower==185);
        }
        if(slowcontrol::Variables::Trigger->HasChanged()) {
            CHECK(firstEvent);
        }
        auto taggerscalers = slowcontrol::Variables::TaggerScalers->Get();
        REQUIRE(taggerscalers.size() == 47);
        if(event.Reconstructed().SlowControls.empty())
            manager.SaveEvent();
        if(!event.Reconstructed().SlowControls.empty())
            CHECK(event.Reconstructed().ID.Lower==395);
        firstEvent = false;
    }

    ~TestPhysics() {
        REQUIRE(nChanged==1);
    }
};



void dotest_ScalerBlobs()
{

    tmpfile_t stage1;
    tmpfile_t stage2;

    // write out some file
    {
        WrapTFileOutput outfile(stage1.filename, true);
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
        WrapTFileOutput outfile(stage2.filename, true);

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

// the processors should behave as follows
//
// S=Skip, B=Buffer, C=Complete, P=Process
//
// TID Timestamp    0 1 2 3 4 5 6 7 8 9 a b c d e f
// TestProcessor1   S S B B B B C B B C B B B B B C
// TestProcessor2   S S S B B B C B B B C B B B B C
// TestProcessor3   C P P P C P P P P P P P P P P C
// TestProcessor4   S S C P P P P P P P P P C P P P

struct procvalue_t {
    procvalue_t(unsigned value, bool hasChanged = false) :
        Value(value), HasChanged(hasChanged)
    {}
    const unsigned Value;
    bool  HasChanged;

    friend ostream& operator<<(ostream& s, const procvalue_t& o) {
        return s << (o.HasChanged ? "'" : " ") << (o.Value==0 ? "S" : to_string(o.Value));
    }

    bool operator==(const procvalue_t& o) const {
        return Value == o.Value && HasChanged == o.HasChanged;
    }

};

struct TestProcessor : slowcontrol::Processor {
    unsigned nProcessed = 0;
    unsigned nPopped = 0;
    unsigned nCompleted = 0;
    queue<unsigned> q;
    virtual return_t ProcessEventData(const TEventData&, physics::manager_t&) override {
        nProcessed++;
        return {};
    }
    virtual void PopQueue() override {
        nPopped++;
        REQUIRE_FALSE(q.empty());
        q.pop();
    }
    void EmplaceQueue() {
        nCompleted++;
        q.emplace(nCompleted);
    }
    virtual unsigned Get() const {
        if(q.empty())
           return 0;
        return q.front();
    }
    virtual vector<procvalue_t> GetExpected() const = 0;

    virtual ~TestProcessor() {
        CHECK(nProcessed==maxEvents);
        CHECK(nCompleted>0);
    }
protected:
    static vector<procvalue_t> MarkChanged(vector<procvalue_t> v) {
        REQUIRE_FALSE(v.empty());
        for(auto it=next(v.begin()); it != v.end(); ++it) {
            if(it->Value != prev(it)->Value)
                it->HasChanged = true;
        }
        return v;
    }
};

struct TestProcessor1 : TestProcessor {

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override
    {
        TestProcessor::ProcessEventData(recon, manager);
        if(recon.ID.Timestamp < 0x2) {
            if(recon.ID.Timestamp == 0x1)
                manager.SaveEvent();
            return return_t::Skip;
        }
        if(recon.ID.Timestamp == 0x6 || recon.ID.Timestamp == 0x9 || recon.ID.Timestamp == 0xf) {
            manager.SaveEvent();
            EmplaceQueue();
            return return_t::Complete;
        }
        return return_t::Buffer;
    }

    virtual vector<procvalue_t> GetExpected() const override {
        return MarkChanged({
         // 0    1    2    3
            {0}, {0}, {1}, {1},
         // 4    5    6    7
            {1}, {1}, {1}, {2},
         // 8    9    a    b
            {2}, {2}, {3}, {3},
         // c    d    e    f
            {3}, {3}, {3}, {3},
        });
    }

    virtual ~TestProcessor1() {
        REQUIRE(nPopped==3);
    }
};

struct TestProcessor2 : TestProcessor {

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override
    {
        TestProcessor::ProcessEventData(recon, manager);
        if(recon.ID.Timestamp < 0x3) {
            if(recon.ID.Timestamp == 0x2)
                manager.SaveEvent();
            return return_t::Skip;
        }
        if(recon.ID.Timestamp == 0x6 || recon.ID.Timestamp == 0xa || recon.ID.Timestamp == 0xf) {
            manager.SaveEvent();
            EmplaceQueue();
            return return_t::Complete;
        }
        return return_t::Buffer;
    }

    virtual vector<procvalue_t> GetExpected() const override {
        return MarkChanged({
         // 0    1    2    3
            {0}, {0}, {0}, {1},
         // 4    5    6    7
            {1}, {1}, {1}, {2},
         // 8    9    a    b
            {2}, {2}, {2}, {3},
         // c    d    e    f
            {3}, {3}, {3}, {3},
        });
    }

    virtual ~TestProcessor2() {
        REQUIRE(nPopped==3);
    }
};

struct TestProcessor3 : TestProcessor {

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) override
    {
        TestProcessor::ProcessEventData(recon, manager);
        if(recon.ID.Timestamp == 0x0 || recon.ID.Timestamp == 0x4 || recon.ID.Timestamp == 0xf) {
            manager.SaveEvent();
            EmplaceQueue();
            return return_t::Complete;
        }
        return return_t::Process;
    }

    virtual vector<procvalue_t> GetExpected() const override {
        return MarkChanged({
         // manually mark the first item as HasChanged
         // 0          1    2    3
            {1, true}, {1}, {1}, {1},
         // 4    5    6    7
            {2}, {2}, {2}, {2},
         // 8    9    a    b
            {2}, {2}, {2}, {2},
         // c    d    e    f
            {2}, {2}, {2}, {3},
        });
    }
    virtual ~TestProcessor3() {
        REQUIRE(nPopped==2);
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
        else if(recon.ID.Timestamp == 0x2 || recon.ID.Timestamp == 0xc) {
            EmplaceQueue();
            return return_t::Complete;
        }
        return return_t::Process;
    }

    virtual vector<procvalue_t> GetExpected() const override {
        return MarkChanged({
         // 0    1    2    3
            {0}, {0}, {1}, {1},
         // 4    5    6    7
            {1}, {1}, {1}, {1},
         // 8    9    a    b
            {1}, {1}, {1}, {1},
         // c    d    e    f
            {2}, {2}, {2}, {2},
        });
    }
    virtual ~TestProcessor4() {
        REQUIRE(nPopped==1);
    }
};

struct TestSlowControlManager : SlowControlManager {
    TestSlowControlManager(const vector<unsigned>& enabled) : SlowControlManager() {
        // previous tests might have requested static slowcontrol variables
        // and the default ctor searches for it...
        processors.clear();
        // we add our own test processors
        if(std_ext::contains(enabled, 1))
            AddProcessor(make_shared<TestProcessor1>());
        if(std_ext::contains(enabled, 2))
            AddProcessor(make_shared<TestProcessor2>());
        if(std_ext::contains(enabled, 3))
            AddProcessor(make_shared<TestProcessor3>());
        if(std_ext::contains(enabled, 4))
            AddProcessor(make_shared<TestProcessor4>());
        CHECK(processors.size() == enabled.size());
    }

    std::vector<std::shared_ptr<TestProcessor>> GetTestProcessors() const {
        std::vector<std::shared_ptr<TestProcessor>> testprocs;
        for(auto& p : processors)
            testprocs.emplace_back(dynamic_pointer_cast<TestProcessor, slowcontrol::Processor>(p.Processor));
        return testprocs;
    }
};

result_t run_TestSlowControlManager(const vector<unsigned>& enabled) {
    TestSlowControlManager scm(enabled);

    // this is basically how PhysicsManager drives the SlowControlManager

    result_t r;

    struct value_t {
        explicit value_t(TID id) : ID(id) {}
        TID ID;
        vector<procvalue_t> ProcValues;

        ostream& operator<<(ostream& s) const {
            return s << hex << ID.Timestamp << dec << " " << ProcValues;
        }

        bool WantsSkip() const {
            for(auto& v : ProcValues)
                if(v.Value == 0) // Value 0 means skip
                    return true;
            return false;
        }

        bool operator==(const value_t& o) const {
            if(ID != o.ID
               || ProcValues.size() != o.ProcValues.size()
               || WantsSkip() != o.WantsSkip())
            {
                return false;
            }
            if(WantsSkip())
                return true;
            // no skip wanted, then ProcValues should match
            return ProcValues == o.ProcValues;
        }
    };

    vector<value_t> values;
    vector<value_t> values_expected;
    while(r.nEventsRead<maxEvents) {
        while(r.nEventsRead<maxEvents) {
            TID tid(r.nEventsRead);
            ++r.nEventsRead;

            input::event_t event;
            event.MakeReconstructed(tid);
            if(scm.ProcessEvent(move(event)))
                break; // became complete, so start popping events
        }

        r.nContextSwitched++;

        while(auto event = scm.PopEvent()) {

            REQUIRE(event.Event.HasReconstructed());
            r.nEventsPopped++;
            r.nEventsSkipped += event.WantsSkip;
            r.nEventsSavedForSC += event.Event.SavedForSlowControls;
            const auto& tid = event.Event.Reconstructed().ID;
            values.emplace_back(tid);
            values_expected.emplace_back(tid);

            for(auto& p : scm.GetTestProcessors()) {
                values.back().ProcValues.emplace_back(event.WantsSkip ? 0 : p->Get(), p->HasChanged());
                values_expected.back().ProcValues.emplace_back( p->GetExpected().at(tid.Timestamp) );
            }
        }
    }

    // propagate hasChanged flag forward to non-skipped events
    // this depends on how many processors are actually activated
    // so it must be done after running the manager
    if(values_expected.size()>1) {
        for(auto it = values_expected.begin();
            it != std::prev(values_expected.end());
            ++it)
        {
            if(it->WantsSkip()) {
                auto& ProcValues = it->ProcValues;
                auto& next_ProcValues = std::next(it)->ProcValues;
                REQUIRE(ProcValues.size() == next_ProcValues.size());
                for(unsigned i=0;i<ProcValues.size();i++) {
                    if(ProcValues[i].HasChanged)
                        next_ProcValues[i].HasChanged = ProcValues[i].HasChanged;
                }
            }
        }
    }

    CHECK(values.size() == values_expected.size());
    CHECK_FALSE(values.empty());
    for(unsigned i=0;i<values.size();i++) {
        CHECK(values[i] == values_expected[i]);
    }

    return r;
}
