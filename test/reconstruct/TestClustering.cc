#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "reconstruct/Clustering.h"

#include "expconfig/ExpConfig.h"

#include "expconfig/detectors/CB.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;


void dotest();

TEST_CASE("Clustering", "[reconstruct]") {
    test::EnsureSetup();
    dotest();
}

void dotest() {
    auto cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    REQUIRE(cb_detector != nullptr);

    // build some readhits
    auto make_readhit = [] (unsigned channel, double energy) {

    };
}
