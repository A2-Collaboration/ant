#include "catch.hpp"
#include "catch_config.h"

#include "reconstruct/Reconstruct.h"
#include "reconstruct/CandidateBuilder.h"
#include "reconstruct/Clustering.h"
#include "reconstruct/UpdateableManager.h"

#include "unpacker/Unpacker.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;


void dotest();

TEST_CASE("CandidateBuilder", "[reconstruct]") {
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

struct counts_t {
    int clusters = 0;
    int candidates = 0;
    int candidateclusters = 0;
    int allclusters = 0;

    counts_t operator-(const counts_t& other) {
        clusters -= other.clusters;
        candidates -= other.candidates;
        candidateclusters -= other.candidateclusters;
        allclusters -= other.allclusters;
        return *this;
    }

};

counts_t getCounts(
        const CandidateBuilder::sorted_clusters_t& sorted_clusters,
        const CandidateBuilder::candidates_t& candidates,
        const CandidateBuilder::clusters_t& all_clusters) {
    counts_t counts;
    counts.clusters = getTotalCount(sorted_clusters);
    counts.candidates = candidates.size();
    counts.allclusters = all_clusters.size();
    for(const auto& cand : candidates)
        counts.candidateclusters += cand->Clusters.size();
    return counts;
}

// we use the friend class trick to test private methods
namespace ant {

struct CandidateBuilderTester : CandidateBuilder {

    using CandidateBuilder::CandidateBuilder; // use base class constructors

    virtual void BuildCandidates(sorted_clusters_t& sorted_clusters,
            candidates_t& candidates,
            clusters_t& all_clusters
            ) override {

        counts_t before;
        counts_t after;
        counts_t diff;

        before = getCounts(sorted_clusters, candidates, all_clusters);
        if(cb && pid)
            Build_PID_CB(sorted_clusters, candidates, all_clusters);
        after = getCounts(sorted_clusters, candidates, all_clusters);
        diff = after - before;
        REQUIRE(diff.candidateclusters + diff.clusters == 0);

        before = getCounts(sorted_clusters, candidates, all_clusters);
        if(taps && tapsveto)
            Build_TAPS_Veto(sorted_clusters, candidates, all_clusters);
        after = getCounts(sorted_clusters, candidates, all_clusters);
        diff = after - before;
        REQUIRE(diff.candidateclusters + diff.clusters == 0);


        before = getCounts(sorted_clusters, candidates, all_clusters);
        Catchall(sorted_clusters, candidates, all_clusters);
        after = getCounts(sorted_clusters, candidates, all_clusters);
        diff = after - before;
        REQUIRE(diff.candidateclusters + diff.clusters == 0);
    }

    virtual void Build(sorted_clusters_t sorted_clusters,
            candidates_t& candidates,
            clusters_t& all_clusters
            ) override {
        counts_t before;
        counts_t after;

        before = getCounts(sorted_clusters, candidates, all_clusters);
        REQUIRE(before.allclusters==0);
        REQUIRE(before.candidateclusters==0);
        REQUIRE(before.candidates==0);
        REQUIRE(before.clusters>0);
        CandidateBuilder::Build(sorted_clusters, candidates, all_clusters);
        after = getCounts(sorted_clusters, candidates, all_clusters);
        REQUIRE(before.clusters == after.allclusters);

        // examine unmatched flag clusters
        size_t unmatched_clusters = 0;
        for(auto cluster : all_clusters) {
            if(cluster->HasFlag(TCluster::Flags_t::Unmatched))
                unmatched_clusters++;
        }
        REQUIRE(unmatched_clusters>0);
        REQUIRE(all_clusters.size() - unmatched_clusters == (unsigned)after.candidateclusters);
    }
};

struct ReconstructTester : Reconstruct {

    virtual void Initialize(const TID& tid) override
    {
        Reconstruct::Initialize(tid);
        // replace the candidate builder with our tester
        const auto& config = ExpConfig::Reconstruct::Get(tid);
        candidatebuilder = std_ext::make_unique<CandidateBuilderTester>(sorted_detectors, config);
    }
};
}


void dotest() {
    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // instead of the usual reconstruct, we use our tester
    ReconstructTester reconstruct;

    unsigned nReads = 0;

    while(auto event = unpacker->NextEvent()) {

        auto& hits = event->Reconstructed->DetectorReadHits;

        if(!hits.empty()) {
            nReads++;
            if(nReads == 6 || nReads == 2) {
                reconstruct.DoReconstruct(*event->Reconstructed);
            }
            else if(nReads > 6)
                break;
        }
    }
}
