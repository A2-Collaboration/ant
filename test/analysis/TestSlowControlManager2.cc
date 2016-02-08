#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "analysis/input/DataReader.h"
#include "analysis/physics/PhysicsManager.h"
#include "analysis/slowcontrol/SlowControlVariables.h"
#include "tree/TSlowControl.h"
#include <chrono>
#include "base/Logger.h"

#include <iostream>
#include <list>

using namespace std;
using namespace ant;
using namespace ant::analysis;

//**** Fake Reader ********************************************************

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
        event.Reconstructed->SlowControls.push_back(s);
    }

    ++i;

    return true;
}

//**** Tools ********************************************************

TSlowControl makeAcqu(const std::string& name, const int value) {
    TSlowControl s(TSlowControl::Type_t::AcquScaler,TSlowControl::Validity_t::Backward,0, name, "Auto Gen");
    s.Payload_Int.push_back(TKeyValue<int64_t>(0,value));
    return s;
}

//**** Tests ********************************************************

void test_FakeReader();

TEST_CASE("SlowControlManager2", "[analysis]") {
    test_FakeReader();
}

void test_FakeReader() {
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

        auto e = TEvent::MakeReconstructed(tid);
        ++tid;

        if(reader.ReadNextEvent(*e)) {
            cout << *e << endl;

        } else {
            break;
        }
    } while(true);

}
