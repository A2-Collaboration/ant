#include "catch.hpp"

#include "base/Detector_t.h"

#include <iostream>

using namespace ant;
using namespace std;

void dotest();

TEST_CASE("Detector_t", "[base]") {
    dotest();
}

const Detector_t::Any_t CBTAPS = Detector_t::Type_t::CB | Detector_t::Type_t::TAPS;

Detector_t::Any_t cb_taps_test(const Detector_t::Any_t& a, const Detector_t::Any_t& b) {
    return (a & CBTAPS) ^ (b & CBTAPS);
}

void dotest() {

    // test printable
    stringstream ss;
    ss << Detector_t::Any_t::CB_Apparatus;
    REQUIRE(ss.str() == "CB|PID|MWPC0|MWPC1");

    // test easy usage in if-statements
    if(Detector_t::Any_t::CB_Apparatus & Detector_t::Type_t::CB) {
        REQUIRE(true);
    }
    else {
        REQUIRE(false);
    }

    REQUIRE((Detector_t::Any_t::CB_Apparatus & Detector_t::Type_t::CB));
    REQUIRE((Detector_t::Type_t::MWPC0 & Detector_t::Any_t::CB_Apparatus));
    REQUIRE((Detector_t::Type_t::MWPC1 & Detector_t::Any_t::CB_Apparatus));
    REQUIRE((Detector_t::Type_t::PID & Detector_t::Any_t::CB_Apparatus));
    REQUIRE_FALSE((Detector_t::Type_t::TAPS & Detector_t::Any_t::CB_Apparatus));
    REQUIRE_FALSE((Detector_t::Type_t::TAPSVeto & Detector_t::Any_t::CB_Apparatus));

    REQUIRE((Detector_t::Type_t::TAPS & Detector_t::Any_t::TAPS_Apparatus));
    REQUIRE((Detector_t::Type_t::TAPSVeto & Detector_t::Any_t::TAPS_Apparatus));

    REQUIRE((Detector_t::Type_t::PID & Detector_t::Any_t::Veto));
    REQUIRE((Detector_t::Type_t::TAPSVeto & Detector_t::Any_t::Veto));

    REQUIRE(Detector_t::Any_t::CB_Apparatus == Detector_t::Any_t::CB_Apparatus);
    REQUIRE(Detector_t::Any_t::CB_Apparatus != Detector_t::Any_t::TAPS_Apparatus);

    Detector_t::Any_t detector1 = Detector_t::Any_t::None;
    detector1 |= Detector_t::Type_t::CB;

    REQUIRE((detector1 & Detector_t::Any_t::CB_Apparatus));
    REQUIRE(detector1 == Detector_t::Type_t::CB);
    REQUIRE(detector1 != Detector_t::Type_t::TAPS);
    REQUIRE(detector1 != Detector_t::Type_t::MWPC0);

    REQUIRE((cb_taps_test(Detector_t::Type_t::CB, Detector_t::Type_t::TAPS) & CBTAPS));
    REQUIRE((cb_taps_test(Detector_t::Type_t::TAPS, Detector_t::Type_t::CB) & CBTAPS));
    REQUIRE_FALSE((cb_taps_test(Detector_t::Type_t::CB, Detector_t::Type_t::CB) & CBTAPS));
    REQUIRE_FALSE((cb_taps_test(Detector_t::Type_t::TAPS, Detector_t::Type_t::TAPS) & CBTAPS));

    REQUIRE((cb_taps_test(Detector_t::Type_t::CB | Detector_t::Type_t::PID, Detector_t::Type_t::TAPS) & CBTAPS));

    REQUIRE(Detector_t::Any_t::CB_Apparatus.test(Detector_t::Type_t::CB));

    Detector_t::ElementFlags_t flags;
    flags.set(Detector_t::ElementFlag_t::BadTDC);
    flags.set(Detector_t::ElementFlag_t::Missing);
    REQUIRE(flags.test(Detector_t::ElementFlag_t::Missing));
    REQUIRE(flags.test(Detector_t::ElementFlag_t::BadTDC));
    flags.unset(Detector_t::ElementFlag_t::Missing);
    REQUIRE_FALSE(flags.test(Detector_t::ElementFlag_t::Missing));
    REQUIRE(flags.test(Detector_t::ElementFlag_t::BadTDC));

    REQUIRE(flags.any());
    REQUIRE_FALSE(flags.none());
    REQUIRE_FALSE(flags.all());
    REQUIRE(static_cast<bool>(flags));

}
