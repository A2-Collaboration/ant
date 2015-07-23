#include "catch.hpp"
#include "catch_config.h"

#include "reconstruct/Reconstruct.h"
#include "reconstruct/TrackBuilder.h"

#include "unpacker/Unpacker.h"

#include "tree/THeaderInfo.h"

#include "base/std_ext.h"

using namespace std;

void dotest();

TEST_CASE("Reconstruct", "[reconstruct]") {
    dotest();
}

// we use the friend class trick to test private methods
namespace ant {
struct ReconstructTester {

    ReconstructTester(const THeaderInfo& headerInfo) :
        r(headerInfo)
    {}

    Reconstruct r;

    unique_ptr<TEvent> DoReconstruct(TDetectorRead& detectorRead)
    {

        // categorize the hits by detector type
        // this is handy for all subsequent reconstruction steps
        map<Detector_t::Type_t, list< TDetectorReadHit* > > sorted_readhits;
        for(TDetectorReadHit& readhit : detectorRead.Hits) {
            sorted_readhits[readhit.GetDetectorType()].push_back(addressof(readhit));
        }

        // apply calibration (this may change the given detectorRead!)
        for(const auto& calib : r.calibrations) {
            calib->ApplyTo(sorted_readhits);
        }

        // for debug purposes, dump out the detectorRead
        //cout << detectorRead << endl;

        // already create the event here, since Tagger
        // doesn't need hit matching and thus can be filled already
        // in BuildHits (see below)
        auto event = std_ext::make_unique<TEvent>(detectorRead.ID);

        // the detectorRead is now calibrated as far as possible
        // lets start the hit matching, which builds the TClusterHit's
        // we also extract the energy, which is always defined as a
        // single value with type Channel_t::Type_t
        Reconstruct::sorted_bydetectortype_t<Reconstruct::HitWithEnergy_t> sorted_clusterhits;
        r.BuildHits(move(sorted_readhits), sorted_clusterhits, event->Tagger);

        // then build clusters (at least for calorimeters this is not trivial)
        Reconstruct::sorted_bydetectortype_t<TCluster> sorted_clusters;
        r.BuildClusters(move(sorted_clusterhits), sorted_clusters);


        // finally, do the track building
        r.trackbuilder->Build(move(sorted_clusters), event);

        // uncomment for debug purposes
        //cout << *event << endl;

        return event;
    }
};
}

using namespace ant;

void dotest() {
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // instead of the usual reconstruct, we use our tester
    unique_ptr<ReconstructTester> reconstruct;

    while(auto item = unpacker->NextItem()) {

        auto HeaderInfo = dynamic_cast<THeaderInfo*>(item.get());
        if(HeaderInfo != nullptr) {
            reconstruct = std_ext::make_unique<ReconstructTester>(*HeaderInfo);
            continue;
        }

        auto DetectorRead = dynamic_cast<TDetectorRead*>(item.get());
        if(DetectorRead != nullptr) {

            if(reconstruct) {
                reconstruct->DoReconstruct(*DetectorRead);
            }
        }
    }
}
