#include "catch.hpp"
#include "catch_config.h"

#include "unpacker/Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "expconfig/detectors/CB.h"

#include <iostream>

using namespace ant;
using namespace std;

void getdetector();
void getlastfound();

TEST_CASE("ExpConfig GetLastFound", "[expconfig]") {
    ExpConfig::Setup::Cleanup();
    getlastfound();
}

TEST_CASE("ExpConfig GetDetector", "[expconfig]") {
    ExpConfig::Setup::Cleanup();
    getdetector();
}

void getdetector() {
    REQUIRE_NOTHROW(Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz"));

    auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    REQUIRE(tagger != nullptr);
    auto tagger_fromtype = ExpConfig::Setup::GetDetector(Detector_t::Type_t::EPT);
    REQUIRE(tagger == tagger_fromtype);

    auto cb = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    REQUIRE(cb != nullptr);
    REQUIRE(cb->GetNChannels() == 720);

    auto ladder = ExpConfig::Setup::GetDetector(Detector_t::Type_t::Tagger);
    REQUIRE(ladder == nullptr);

    // switch to another setup
    REQUIRE_NOTHROW(Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-small.dat.xz"));
    auto tagger_notthere = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    REQUIRE(tagger_notthere == nullptr);
}

void getlastfound() {
    REQUIRE(ExpConfig::Setup::GetLastFound() == nullptr);
    ExpConfig::Setup::ManualName = "Setup_2014_07_EPT_Prod";
    REQUIRE(ExpConfig::Setup::GetLastFound() != nullptr);
}
