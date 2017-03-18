#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "unpacker/Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "expconfig/detectors/CB.h"

#include <iostream>

using namespace ant;
using namespace std;

void getdetector();
void getlastfound();
void getall();

TEST_CASE("ExpConfig Get (all)", "[expconfig]") {
    getall();
}

TEST_CASE("ExpConfig GetDetector", "[expconfig]") {

    getdetector();
}

void getall() {
    auto setupnames = ExpConfig::Setup::GetNames();
    for(auto setupname : setupnames) {
        REQUIRE_NOTHROW(ExpConfig::Setup::SetByName(setupname));
    }
}

void getdetector() {
    test::EnsureSetup();
    REQUIRE_NOTHROW(Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz"));

    auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    REQUIRE(tagger.get() != nullptr);
    auto tagger_fromtype = ExpConfig::Setup::GetDetector(Detector_t::Type_t::EPT);
    REQUIRE(tagger == tagger_fromtype);

    auto cb = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    REQUIRE(cb.get() != nullptr);
    REQUIRE(cb->GetNChannels() == 720);

    REQUIRE_THROWS_AS(ExpConfig::Setup::GetDetector(Detector_t::Type_t::Tagger), ExpConfig::Exception);
}

