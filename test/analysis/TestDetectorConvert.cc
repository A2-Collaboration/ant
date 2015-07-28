#include "catch.hpp"

#include "analysis/input/ant/detail/Convert.h"
#include "analysis/Detector.h"
#include "expconfig/Detector_t.h"

#include <iostream>

using namespace std;
using namespace ant;

void dotest();
void dotest_patterns();

const auto CBTAPS = detector_t::CB | detector_t::TAPS;

TEST_CASE("DetectorConvert", "[analysis]") {
    dotest();
}

void dotest() {

    REQUIRE((detector_t(Detector_t::ToBitfield(Detector_t::Type_t::CB)) & detector_t::CB));
    REQUIRE((detector_t(Detector_t::ToBitfield(Detector_t::Type_t::TAPS)) & detector_t::TAPS));
    REQUIRE_FALSE((detector_t(Detector_t::ToBitfield(Detector_t::Type_t::TAPS)) & detector_t::CB));

}

TEST_CASE("Patterns","[analysis]") {
    dotest_patterns();
}

detector_t cb_taps_test(const detector_t& a, const detector_t& b) {
    return (a & CBTAPS) ^ (b & CBTAPS);
}

void dotest_patterns()
{
    // req. 1 cb and 1 taps
    REQUIRE((cb_taps_test(detector_t::CB, detector_t::TAPS) & CBTAPS));
    REQUIRE((cb_taps_test(detector_t::TAPS, detector_t::CB) & CBTAPS));
    REQUIRE_FALSE((cb_taps_test(detector_t::CB, detector_t::CB) & CBTAPS));
    REQUIRE_FALSE((cb_taps_test(detector_t::TAPS, detector_t::TAPS) & CBTAPS));

    REQUIRE((cb_taps_test(detector_t::CB|detector_t::PID, detector_t::TAPS) & CBTAPS));

}
