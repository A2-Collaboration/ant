#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "Unpacker.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include <iostream>
#include <string>

using namespace std;
using namespace ant;

void dotest_problematic();
void dotest_another();

TEST_CASE("Test UnpackerAcqu: Mk1 with error scaler block", "[unpacker]") {
    dotest_problematic();
}

TEST_CASE("Test UnpackerAcqu: Mk1 with another scaler block", "[unpacker]") {
    dotest_another();
}

void dotest_problematic() {
    ant::test::EnsureSetup();
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/AcquMk1_problematic.dat.gz");

    unsigned nSlowControls = 0;
    unsigned nEvents = 0;
    unsigned nHits = 0;
    unsigned nEmptyEvents = 0;

    while(auto event = unpacker->NextEvent()) {
        auto& readhits = event.Reconstructed().DetectorReadHits;
        nEvents++;
        nHits += readhits.size();

        // last event should report proper end of file
        if(nEvents==85) {
            REQUIRE(event.Reconstructed().UnpackerMessages.size()==1);
            REQUIRE(event.Reconstructed().UnpackerMessages.back().Message == "Found proper end of file");
        }

        auto& slowcontrols = event.Reconstructed().SlowControls;
        if(!slowcontrols.empty()) {
            nSlowControls += slowcontrols.size();
        }
    }

    CHECK(nSlowControls == 0);
    CHECK(nEvents == 85);
    CHECK(nHits == 5737);
    CHECK(nEmptyEvents == 0);
}

void dotest_another() {
    ant::test::EnsureSetup();
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/AcquMk1_scalerblock.dat.xz");

    unsigned nSlowControls = 0;
    unsigned nEvents = 0;
    unsigned nHits = 0;
    unsigned nEmptyEvents = 0;

    while(auto event = unpacker->NextEvent()) {
        auto& readhits = event.Reconstructed().DetectorReadHits;
        nEvents++;
        nHits += readhits.size();

        // last event should report proper end of file
        if(nEvents==37) {
            REQUIRE(event.Reconstructed().UnpackerMessages.size()==1);
            REQUIRE(event.Reconstructed().UnpackerMessages.back().Message == "Found proper end of file");
        }

        auto& slowcontrols = event.Reconstructed().SlowControls;
        if(!slowcontrols.empty()) {
            nSlowControls += slowcontrols.size();
        }
    }

    CHECK(nSlowControls == 13);
    CHECK(nEvents == 37);
    CHECK(nHits == 1102);
    CHECK(nEmptyEvents == 0);
}