#include "catch.hpp"
#include "catch_config.h"

#include "reconstruct/Reconstruct.h"
#include "reconstruct/TrackBuilder.h"
#include "reconstruct/Clustering.h"

#include "unpacker/Unpacker.h"

#include "tree/THeaderInfo.h"

#include "base/std_ext.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;


void dotest();

TEST_CASE("Reconstruct", "[reconstruct]") {
    dotest();
}

template<typename T>
unsigned getTotalCount(const Reconstruct::sorted_bydetectortype_t<T>& m) {
    unsigned total = 0;
    for(const auto& m_item : m) {
        const list<T>& list = m_item.second;
        total += list.size();
    }
    return total;
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
        /// \todo Improve requirements


        // apply the calibrations,
        CalibrationApply_traits::readhits_t sorted_readhits;
        r.ApplyCalibrations(detectorRead, sorted_readhits);
        size_t n_readhits = getTotalCount(sorted_readhits);
        REQUIRE(n_readhits>0);

        // already create the event here, since Tagger
        // doesn't need hit matching and thus can be filled already
        // in BuildHits (see below)
        auto event = std_ext::make_unique<TEvent>(detectorRead.ID);

        // the detectorRead is now calibrated as far as possible
        // lets start the hit matching, which builds the TClusterHit's
        // we also extract the energy, which is always defined as a
        // single value with type Channel_t::Type_t
        Reconstruct::sorted_bydetectortype_t<AdaptorTClusterHit> sorted_clusterhits;
        r.BuildHits(move(sorted_readhits), sorted_clusterhits, event->Tagger);
        size_t n_clusterhits = getTotalCount(sorted_clusterhits);
        REQUIRE(n_clusterhits + event->Tagger.Hits.size() <= n_readhits);

        // then build clusters (at least for calorimeters this is not trivial)
        Reconstruct::sorted_bydetectortype_t<TCluster> sorted_clusters;
        r.BuildClusters(move(sorted_clusterhits), sorted_clusters, event->InsaneClusters);
        size_t n_clusters = getTotalCount(sorted_clusters);
        REQUIRE(n_clusters>0);
        REQUIRE(n_clusters <= n_clusterhits);


        // finally, do the track building
        r.trackbuilder->Build(move(sorted_clusters), event->Tracks);
        REQUIRE(!event->Tracks.empty());
        REQUIRE(event->Tracks.size() <= n_clusters);

        return event;
    }
};
}


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
