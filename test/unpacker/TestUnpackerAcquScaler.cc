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

void dotest();

TEST_CASE("Test UnpackerAcqu: Scaler block", "[unpacker]") {
    dotest();
}

void dotest() {
    ant::test::EnsureSetup();
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_scalerblock.dat.xz");

    unsigned nSlowControls = 0;
    unsigned nEvents = 0;
    unsigned nHits = 0;
    unsigned nEmptyEvents = 0;

    bool taggerScalerBlockFound = false;

    while(auto event = unpacker->NextEvent()) {
        auto& readhits = event->Reconstructed().DetectorReadHits;
        nEvents++;
        nHits += readhits.size();

        // last event should report proper end of file
        if(nEvents==211) {
            REQUIRE(event->Reconstructed().UnpackerMessages.back().Message == "Found proper end of file");
        }

        auto& slowcontrols = event->Reconstructed().SlowControls;
        if(!slowcontrols.empty()) {
            nSlowControls += slowcontrols.size();
            for(auto& sc : slowcontrols) {
                // the test setup and the Acqu test blob uses the end-point tagger,
                // so there should be corresponding scalers
                if(sc.Name == "EPT_Scalers") {
                    taggerScalerBlockFound = true;
                    // we know the file is extracted from an EPT run...
                    REQUIRE(sc.Payload_Int.size() == 47);
                }
            }
        }
    }

    REQUIRE(nSlowControls == 15);
    REQUIRE(nEvents == 211);
    REQUIRE(nHits == 28665);
    REQUIRE(nEmptyEvents == 0);
    REQUIRE(taggerScalerBlockFound);
}
