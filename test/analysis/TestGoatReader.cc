#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "analysis/input/goat/GoatReader.h"

#include "base/WrapTFile.h"
#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis::input;
void dotest_read();

TEST_CASE("GoatReader: Read some events", "[analysis]") {
    test::EnsureSetup();
    dotest_read();
}

void dotest_read() {
    auto inputfiles = make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/GoAT_5711_100events.root");

    GoatReader reader(inputfiles);

    REQUIRE(reader.IsSource());

    unsigned nEvents = 0;
    unsigned nCandidates = 0;
    while(true) {
        event_t event;

        if(!reader.ReadNextEvent(event))
            break;

        REQUIRE(event.HasReconstructed());
        REQUIRE_FALSE(event.Reconstructed().ID.IsInvalid());

        nEvents++;
        nCandidates += event.Reconstructed().Candidates.size();
    }

    CHECK(nEvents==100);
    CHECK(nCandidates == 345);
}


