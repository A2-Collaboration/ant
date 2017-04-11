#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "reconstruct/Clustering.h"
#include "reconstruct/Reconstruct.h"
#include "unpacker/Unpacker.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "expconfig/ExpConfig.h"

#include "expconfig/detectors/CB.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;


void dotest_build();
void dotest_statistical();

TEST_CASE("Clustering: Build", "[reconstruct]") {
    test::EnsureSetup();
    dotest_build();
}

TEST_CASE("Clustering: Statistical", "[reconstruct]") {
    test::EnsureSetup();
    dotest_statistical();
}


void dotest_build() {
    auto cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    REQUIRE(cb_detector != nullptr);

    /// \todo Implement building detector hits
}

struct ClusteringTester : Clustering_NextGen {

    using Clustering_NextGen::Clustering_NextGen;

    void Build(const ClusterDetector_t& clusterdetector,
               const TClusterHitList& clusterhits,
               TClusterList& clusters
               ) const override
    {


        double total_energy_before = 0.0;
        for(const TClusterHit& clusterhit : clusterhits) {
            if(clusterhit.IsSane())
                total_energy_before += clusterhit.Energy;
        }

        REQUIRE(clusters.empty());

        Clustering_NextGen::Build(clusterdetector, clusterhits, clusters);

        double total_energy_after = 0.0;
        for(const TCluster& cluster : clusters)
            total_energy_after += cluster.Energy;
        REQUIRE(total_energy_after == Approx(total_energy_before));
    }

};

struct ReconstructTester : Reconstruct {
    ReconstructTester() :
        Reconstruct(std_ext::make_unique<ClusteringTester>())
    {}
};

void dotest_statistical() {
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // instead of the usual reconstruct, we use our tester
    ReconstructTester reconstruct;

    unsigned nClusters = 0;
    unsigned nClusterHits = 0;
    unsigned nSplitClusters = 0;
    unsigned nSplitClusterHits = 0;
    unsigned nTouchesHole_CB = 0;
    unsigned nTouchesHole_TAPS = 0;
    unsigned nTouchesHoleCrystal_CB = 0;
    unsigned nTouchesHoleCrystal_TAPS = 0;


    while(auto event = unpacker->NextEvent()) {
        auto& recon = event.Reconstructed();
        reconstruct.DoReconstruct(recon);
        for(const TCluster& cluster : recon.Clusters) {
            nClusters++;
            nClusterHits += cluster.Hits.size();
            if(cluster.HasFlag(TCluster::Flags_t::Split)) {
                nSplitClusters++;
                nSplitClusterHits += cluster.Hits.size();
            }
            if(cluster.HasFlag(TCluster::Flags_t::TouchesHoleCentral)) {
                REQUIRE(cluster.HasFlag(TCluster::Flags_t::TouchesHoleCrystal));
                if(cluster.DetectorType == Detector_t::Type_t::CB)
                    nTouchesHole_CB++;
                else if(cluster.DetectorType == Detector_t::Type_t::TAPS)
                    nTouchesHole_TAPS++;
            }
            if(cluster.HasFlag(TCluster::Flags_t::TouchesHoleCrystal)) {
                if(cluster.DetectorType == Detector_t::Type_t::CB)
                    nTouchesHoleCrystal_CB++;
                else if(cluster.DetectorType == Detector_t::Type_t::TAPS)
                    nTouchesHoleCrystal_TAPS++;
            }
        }
    }

    CHECK(nClusters == 3279);
    CHECK(nClusterHits == 6068);
    CHECK(nSplitClusters == 133);
    CHECK(nSplitClusterHits == 770);
    CHECK(nTouchesHole_CB == 229);
    CHECK(nTouchesHole_TAPS == 96);
    CHECK(nTouchesHoleCrystal_CB == 314);
    CHECK(nTouchesHoleCrystal_TAPS == 99);
}
