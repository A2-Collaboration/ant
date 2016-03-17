#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "reconstruct/Reconstruct.h"
#include "reconstruct/CandidateBuilder.h"
#include "reconstruct/Clustering.h"
#include "reconstruct/UpdateableManager.h"

#include "unpacker/Unpacker.h"


using namespace std;
using namespace ant;
using namespace ant::reconstruct;


void dotest();

TEST_CASE("Reconstruct", "[reconstruct]") {
    test::EnsureSetup();
    dotest();
}

template<typename T>
unsigned getTotalCount(const T& m) {
    unsigned total = 0;
    for(const auto& m_item : m) {
        const auto& list = m_item.second;
        total += list.size();
    }
    return total;
}

struct ReconstructTester : Reconstruct {

    void DoReconstruct(TEventData& reconstructed) override
    {
        // ignore empty events
        if(reconstructed.DetectorReadHits.empty())
            return;

        /// \todo Improve requirements

        if(!initialized) {
            Initialize(reconstructed.ID);
        }
        REQUIRE(initialized);

        // update the updateables :)
        updateablemanager->UpdateParameters(reconstructed.ID);

        // apply the hooks (mostly calibrations)
        ApplyHooksToReadHits(reconstructed.DetectorReadHits);
        // manually scan the r.sorted_readhits
        // they are a member variable for performance reasons
        size_t n_readhits = 0;
        for(const auto& readhit : sorted_readhits ) {
            n_readhits += readhit.second.size();
        }
        if(!reconstructed.DetectorReadHits.empty())
            REQUIRE(n_readhits>0);

        // the detectorRead is now calibrated as far as possible
        // lets start the hit matching, which builds the TClusterHit's
        // we also extract the energy, which is always defined as a
        // single value with type Channel_t::Type_t
        Reconstruct::sorted_bydetectortype_t<TClusterHit> sorted_clusterhits;
        BuildHits(sorted_clusterhits, reconstructed.TaggerHits);

        // apply hooks which modify clusterhits
        for(const auto& hook : hooks_clusterhits) {
            hook->ApplyTo(sorted_clusterhits);
        }

        size_t n_clusterhits = getTotalCount(sorted_clusterhits);
        REQUIRE(n_clusterhits + reconstructed.TaggerHits.size() <= n_readhits);

        // then build clusters (at least for calorimeters this is not trivial)
        Reconstruct::sorted_clusters_t sorted_clusters;
        BuildClusters(move(sorted_clusterhits), sorted_clusters);
        size_t n_clusters = getTotalCount(sorted_clusters);
        if(!reconstructed.DetectorReadHits.empty())
            REQUIRE(n_clusters>0);
        REQUIRE(n_clusters <= n_clusterhits);

        // apply hooks which modify clusters
        for(const auto& hook : hooks_clusters) {
            hook->ApplyTo(sorted_clusters);
        }

        // do the candidate building
        const auto n_all_before = reconstructed.Clusters.size();
        REQUIRE(n_all_before==0);
        candidatebuilder->Build(move(sorted_clusters),
                                reconstructed.Candidates, reconstructed.Clusters);

        // apply hooks which may modify the whole event
        for(const auto& hook : hooks_eventdata) {
            hook->ApplyTo(reconstructed);
        }

        const auto n_all_after = reconstructed.Clusters.size();
        REQUIRE(n_all_after>=n_all_before);
        const auto n_all_added = n_all_after - n_all_before;
        REQUIRE(n_all_added == n_clusters);

        bool matched_clusters = false;
        for(auto& cluster : reconstructed.Clusters) {
            if(!cluster.HasFlag(TCluster::Flags_t::Unmatched)) {
                matched_clusters = true;
                break;
            }
        }
        if(matched_clusters) {
            const size_t n_candidates = reconstructed.Candidates.size();
            REQUIRE(n_candidates>0);
            REQUIRE(reconstructed.Trigger.CBEnergySum > 0);
        }
    }
};

void dotest() {
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // instead of the usual reconstruct, we use our tester
    ReconstructTester reconstruct;

    unsigned nReads = 0;
    unsigned nHits = 0;
    unsigned nCandidates = 0;

    while(auto event = unpacker->NextEvent()) {

        auto& hits = event.Reconstructed().DetectorReadHits;
        nHits += hits.size();
        if(!hits.empty())
            nReads++;
        reconstruct.DoReconstruct(event.Reconstructed());
        nCandidates += event.Reconstructed().Candidates.size();

    }

    REQUIRE(nReads == 221);
    REQUIRE(nHits == 30260);
    REQUIRE(nCandidates == 860);


}
