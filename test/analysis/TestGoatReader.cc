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
    unsigned nClusters = 0;
    unsigned nCandidates = 0;
    unsigned nTaggerHits = 0;
    double SumCBESum = 0;
    while(true) {
        event_t event;

        if(!reader.ReadNextEvent(event))
            break;

        REQUIRE(event.HasReconstructed());
        auto& recon = event.Reconstructed();
        REQUIRE_FALSE(recon.ID.IsInvalid());
        REQUIRE(recon.ID.isSet(TID::Flags_t::AdHoc));

        nEvents++;
        nCandidates += recon.Candidates.size();
        nClusters += recon.Clusters.size();
        nTaggerHits += recon.TaggerHits.size();

        SumCBESum += recon.Trigger.CBEnergySum;
    }

    CHECK(nEvents==100);
    CHECK(nCandidates == 345);
    CHECK(nClusters == 518);
    CHECK(nTaggerHits == 2831);
    CHECK(SumCBESum/nEvents == Approx(720.344472).epsilon(0.0001));
}


