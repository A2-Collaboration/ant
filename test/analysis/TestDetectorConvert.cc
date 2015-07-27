#include "catch.hpp"

#include "analysis/input/ant/detail/Convert.h"
#include "analysis/Detector.h"
#include "expconfig/Detector_t.h"

#include <iostream>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("DetectorConvert", "[analysis]") {
    dotest();
}

void dotest() {

    REQUIRE((detector_t(Detector_t::ToBitfield(Detector_t::Type_t::CB)) & detector_t::CB));
    REQUIRE((detector_t(Detector_t::ToBitfield(Detector_t::Type_t::TAPS)) & detector_t::TAPS));
    REQUIRE_FALSE((detector_t(Detector_t::ToBitfield(Detector_t::Type_t::TAPS)) & detector_t::CB));

}
