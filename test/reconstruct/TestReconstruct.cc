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
struct ReconstructTester : Reconstruct_traits {
    Reconstruct r;

    void Initialize(const THeaderInfo& headerInfo) override
    {
        r.Initialize(headerInfo);
    }

    MemoryPool<TEvent>::Item DoReconstruct(TDetectorRead& detectorRead) override
    {
        /// \todo Improve requirements

        // update the updateables :)
        r.updateablemanager->UpdateParameters(detectorRead.ID);

        // apply the hooks (mostly calibrations)
        r.ApplyHooksToReadHits(detectorRead);
        // manually scan the r.sorted_readhits
        // they are a member variable for performance reasons
        size_t n_readhits = 0;
        for(const auto& readhit : r.sorted_readhits ) {
            n_readhits += readhit.second.size();
        }
        REQUIRE(n_readhits>0);

        // already create the event here, since Tagger
        // doesn't need hit matching and thus can be filled already
        // in BuildHits (see below)
        auto event = MemoryPool<TEvent>::Get();
        event->ID = detectorRead.ID;

        // the detectorRead is now calibrated as far as possible
        // lets start the hit matching, which builds the TClusterHit's
        // we also extract the energy, which is always defined as a
        // single value with type Channel_t::Type_t
        Reconstruct::sorted_bydetectortype_t<AdaptorTClusterHit> sorted_clusterhits;
        r.BuildHits(sorted_clusterhits, event->Tagger);

        // apply hooks which modify clusterhits
        for(const auto& hook : r.hooks_clusterhits) {
            hook->ApplyTo(sorted_clusterhits);
        }

        size_t n_clusterhits = getTotalCount(sorted_clusterhits);
        REQUIRE(n_clusterhits + event->Tagger.Hits.size() <= n_readhits);

        // then build clusters (at least for calorimeters this is not trivial)
        Reconstruct::sorted_bydetectortype_t<TCluster> sorted_clusters;
        r.BuildClusters(move(sorted_clusterhits), sorted_clusters);
        size_t n_clusters = getTotalCount(sorted_clusters);
        REQUIRE(n_clusters>0);
        REQUIRE(n_clusters <= n_clusterhits);

        // apply hooks which modify clusters
        for(const auto& hook : r.hooks_clusters) {
            hook->ApplyTo(sorted_clusters);
        }

        // finally, do the candidate building
        const auto n_all_before = event->AllClusters.size();
        REQUIRE(n_all_before==0);
        r.candidatebuilder->Build(move(sorted_clusters), event->Candidates, event->AllClusters);

        const auto n_all_after = event->AllClusters.size();
        REQUIRE(n_all_after>=n_all_before);
        const auto n_all_added = n_all_after - n_all_before;
        REQUIRE(n_all_added == n_clusters);

        bool matched_clusters = false;
        for(auto cluster : event->AllClusters) {
            if(!cluster.HasFlag(TCluster::Flags_t::Unmatched)) {
                matched_clusters = true;
                break;
            }
        }
        if(matched_clusters) {
            const size_t n_candidates = event->Candidates.size();
            REQUIRE(n_candidates>0);
        }

        return event;
    }
};
}


void dotest() {
    ExpConfig::Setup::ManualName = "Setup_Test";

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
            if(reconstruct) {
                auto event = reconstruct->DoReconstruct(*DetectorRead);
                nCandidates += event->Candidates.size();
            }
        }
    }

    REQUIRE(nReads == 221);
    REQUIRE(nHits == 30039);
    REQUIRE(nCandidates == 822);


}
