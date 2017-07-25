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

    REQUIRE((reader.GetFlags() & reader_flag_t::IsSource));

    unsigned nEvents = 0;
    unsigned nClusters = 0;
    unsigned nCandidates = 0;
    unsigned nTaggerHits = 0;
    double SumCBESum = 0;
    double SumMultiplicity = 0;
    double SumCBTiming = 0;

    map<Detector_t::Type_t, unsigned> readHitsByDetector;
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
        SumMultiplicity += recon.Trigger.ClusterMultiplicity;

        CHECK(recon.Trigger.DAQEventID != 0);
        CHECK(recon.Trigger.DAQErrors.empty()); // no errors read-out

        for(auto& readhit : recon.DetectorReadHits) {
            readHitsByDetector[readhit.DetectorType]++;
        }
    }

    CHECK(nEvents==100);
    CHECK(nCandidates == 345);
    CHECK(nClusters == 518);
    CHECK(nTaggerHits == 2831);
    CHECK(SumCBESum/nEvents == Approx(720.344472).epsilon(0.0001));
    // multiplicity is always zero in input blob....
    CHECK(SumMultiplicity/nEvents == Approx(0).epsilon(0.0001));

    // no EPT/Tagger readhits given by Acqu/TA2GoAT
    CHECK(readHitsByDetector.size() == 4);
    CHECK(readHitsByDetector[Detector_t::Type_t::CB] == 2701);
    CHECK(readHitsByDetector[Detector_t::Type_t::TAPS] == 154);
    CHECK(readHitsByDetector[Detector_t::Type_t::PID] == 622);
    CHECK(readHitsByDetector[Detector_t::Type_t::TAPSVeto] == 136);
}


