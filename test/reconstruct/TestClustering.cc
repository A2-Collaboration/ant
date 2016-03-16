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

    // build some readhits
//    auto make_readhit = [] (unsigned channel, double energy) {

//    };
}

struct ClusteringTester : Clustering {

    using Clustering::Clustering;

    void Build(const std::shared_ptr<ClusterDetector_t>& clusterdetector,
               const TClusterHitList& clusterhits,
               TClusterList& clusters
               ) override
    {


        double total_energy_before = 0.0;
        for(const TClusterHit& clusterhit : clusterhits) {
            if(clusterhit.IsSane())
                total_energy_before += clusterhit.Energy;
        }

        REQUIRE(clusters.empty());

        Clustering::Build(clusterdetector, clusterhits, clusters);

        double total_energy_after = 0.0;
        for(const TCluster& cluster : clusters)
            total_energy_after += cluster.Energy;
        REQUIRE(total_energy_after == Approx(total_energy_before));
    }

};

struct ReconstructTester : Reconstruct {

    void Initialize(const TID& tid) override {
        Reconstruct::Initialize(tid);
        // replace the clustering with our tester
        const auto& config = ExpConfig::Reconstruct::Get(tid);
        clustering = std_ext::make_unique<ClusteringTester>(config);
    }

};

void dotest_statistical() {
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // instead of the usual reconstruct, we use our tester
    ReconstructTester reconstruct;

    unsigned nClusters = 0;
    unsigned nClusterHits = 0;
    unsigned nSplitClusters = 0;
    unsigned nSplitClusterHits = 0;

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

        }
    }

    REQUIRE(nClusters == 3251);
    REQUIRE(nClusterHits == 5792);
    REQUIRE(nSplitClusters == 107);
    REQUIRE(nSplitClusterHits == 671);
}
