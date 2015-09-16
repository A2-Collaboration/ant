#include "catch.hpp"
#include "catch_config.h"

#include "reconstruct/Reconstruct.h"
#include "reconstruct/CandidateBuilder.h"
#include "reconstruct/Clustering.h"
#include "reconstruct/UpdateableManager.h"

#include "unpacker/Unpacker.h"

#include "tree/THeaderInfo.h"
#include "tree/TDetectorRead.h"

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
    int insaneclusters = 0;

    counts_t operator-(const counts_t& other) {
        clusters -= other.clusters;
        candidates -= other.candidates;
        candidateclusters -= other.candidateclusters;
        insaneclusters -= other.insaneclusters;
        return *this;
    }

};

counts_t getCounts(
        const std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters,
        const TEvent::candidates_t& candidates,
        const std::vector<TCluster>& insane_clusters) {
    counts_t counts;
    counts.clusters = getTotalCount(sorted_clusters);
    counts.candidates = candidates.size();
    counts.insaneclusters = insane_clusters.size();
    for(const auto& cand : candidates)
        counts.candidateclusters += cand.Clusters.size();
    return counts;
}

// we use the friend class trick to test private methods
namespace ant {

struct CandidateBuilderTester : CandidateBuilder {

    using CandidateBuilder::CandidateBuilder; // use base class constructors

    virtual void Build(std::map<Detector_t::Type_t, std::list<TCluster> > sorted_clusters,
            TEvent::candidates_t& candidates,
            std::vector<TCluster>& insane_clusters
            ) {

        counts_t before;
        counts_t after;
        counts_t diff;

        before = getCounts(sorted_clusters, candidates, insane_clusters);
        if(cb && pid)
            Build_PID_CB(sorted_clusters, candidates);
        after = getCounts(sorted_clusters, candidates, insane_clusters);
        diff = after - before;
        REQUIRE(diff.candidateclusters + diff.clusters == 0);

        before = getCounts(sorted_clusters, candidates, insane_clusters);
        if(taps && tapsveto)
            Build_TAPS_Veto(sorted_clusters, candidates);
        after = getCounts(sorted_clusters, candidates, insane_clusters);
        diff = after - before;
        REQUIRE(diff.candidateclusters + diff.clusters == 0);


        before = getCounts(sorted_clusters, candidates, insane_clusters);

        Catchall(sorted_clusters, candidates);

        // move the rest to insane clusters
        for(auto& det_entry : sorted_clusters) {
            for(auto& cluster : det_entry.second) {
                insane_clusters.emplace_back(cluster);
            }
        }
        sorted_clusters.clear();

        after = getCounts(sorted_clusters, candidates, insane_clusters);
        diff = after - before;

        REQUIRE(diff.candidateclusters == diff.candidates);
        REQUIRE(diff.insaneclusters + diff.candidateclusters + diff.clusters == 0);

    }
};

struct ReconstructTester : Reconstruct {

    virtual void Initialize(const THeaderInfo& headerInfo) override
    {
        Reconstruct::Initialize(headerInfo);
        // replace the candidate builder with our tester
        const auto& config = ExpConfig::Reconstruct::Get(headerInfo);
        candidatebuilder = std_ext::make_unique<CandidateBuilderTester>(sorted_detectors, config);
    }
};
}


void dotest() {
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // instead of the usual reconstruct, we use our tester
    unique_ptr<ReconstructTester> reconstruct;

    unsigned nReads = 0;
    unsigned nHits = 0;
    unsigned nCandidates = 0;

    while(auto item = unpacker->NextItem()) {

        auto HeaderInfo = dynamic_cast<THeaderInfo*>(item.get());
        if(HeaderInfo != nullptr) {
            reconstruct = std_ext::make_unique<ReconstructTester>();
            reconstruct->Initialize(*HeaderInfo);
            continue;
        }

        auto DetectorRead = dynamic_cast<TDetectorRead*>(item.get());
        if(DetectorRead != nullptr) {
            nReads++;
            nHits += DetectorRead->Hits.size();
            if(reconstruct && (nReads == 6 || nReads == 2)) {
                auto event = reconstruct->DoReconstruct(*DetectorRead);
                nCandidates += event->Candidates.size();
            }
            else if(nReads > 6)
                break;
        }
    }
}
