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


void dotest_sanity();
void dotest_ignoredelements_raw();
void dotest_ignoredelements_raw_include();
void dotest_ignoredelements_geant();
void dotest_ignoredelements_geant_include();


TEST_CASE("Reconstruct: Chain sanity checks", "[reconstruct]") {
    test::EnsureSetup();
    dotest_sanity();
}

TEST_CASE("Reconstruct: Ignored elements raw data", "[reconstruct]") {
    test::EnsureSetup();
    dotest_ignoredelements_raw();
}

TEST_CASE("Reconstruct: Include Ignored elements with raw data", "[reconstruct]") {
    test::EnsureSetup(true);
    dotest_ignoredelements_raw_include();
}

TEST_CASE("Reconstruct: Ignored elements geant", "[reconstruct]") {
    test::EnsureSetup();
    dotest_ignoredelements_geant();
}

TEST_CASE("Reconstruct: Include ignored elements with geant", "[reconstruct]") {
    test::EnsureSetup(true);
    dotest_ignoredelements_geant_include();
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

    void DoReconstruct(TEventData& reconstructed) const override
    {
        // ignore empty events
        if(reconstructed.DetectorReadHits.empty())
            return;

        /// \todo Improve requirements

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
        }
    }
};

void dotest_sanity() {
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // instead of the usual reconstruct, we use our tester
    ReconstructTester reconstruct;

    unsigned nReads = 0;
    unsigned nHits = 0;
    unsigned nCandidates = 0;
    unsigned nCandidatesCBPID = 0;
    unsigned nCandidatesTAPSVeto = 0;


    while(auto event = unpacker->NextEvent()) {

        auto& hits = event.Reconstructed().DetectorReadHits;
        nHits += hits.size();
        if(!hits.empty())
            nReads++;
        reconstruct.DoReconstruct(event.Reconstructed());
        nCandidates += event.Reconstructed().Candidates.size();
        for(auto& cand : event.Reconstructed().Candidates) {
            if(cand.Detector & Detector_t::Type_t::CB &&
               cand.Detector & Detector_t::Type_t::PID)
                nCandidatesCBPID++;
            if(cand.Detector & Detector_t::Type_t::TAPS &&
               cand.Detector & Detector_t::Type_t::TAPSVeto)
                nCandidatesTAPSVeto++;

        }
    }

    CHECK(nReads == 221);
    CHECK(nHits == 32243);
    CHECK(nCandidates == 864);
    CHECK(nCandidatesCBPID == 276);
    CHECK(nCandidatesTAPSVeto == 146);

}

map<Detector_t::Type_t, unsigned> getReconstructedHits(bool geant) {
    auto unpacker = Unpacker::Get(geant ?  string(TEST_BLOBS_DIRECTORY)+"/Geant_with_TID.root" :
                                           string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
    Reconstruct reconstruct;
    map<Detector_t::Type_t, unsigned> hits;
    auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    while(auto event = unpacker->NextEvent()) {
        reconstruct.DoReconstruct(event.Reconstructed());
        for(auto& cluster : event.Reconstructed().Clusters) {
            hits[cluster.DetectorType] += cluster.Hits.size();
        }
        for(auto& taggerhit : event.Reconstructed().TaggerHits) {
            hits[tagger->Type] += taggerhit.Electrons.size();
        }
    }
    return hits;
}

void dotest_ignoredelements_raw() {
    auto clusterHits_before = getReconstructedHits(false);
    CHECK(clusterHits_before.size() == 5);
    CHECK(clusterHits_before[Detector_t::Type_t::CB] == 3759);
    CHECK(clusterHits_before[Detector_t::Type_t::TAPS] == 720);
    CHECK(clusterHits_before[Detector_t::Type_t::PID] == 633);
    CHECK(clusterHits_before[Detector_t::Type_t::TAPSVeto] == 956);
    CHECK(clusterHits_before[Detector_t::Type_t::EPT] == 6272);

    // ignore every second element in each detector type
    for(const auto& item : clusterHits_before) {
        auto det = ExpConfig::Setup::GetDetector(item.first);
        for(unsigned ch=0;ch<det->GetNChannels();ch++) {
            if(ch % 2 == 0)
                det->SetElementFlags(ch, Detector_t::ElementFlag_t::Broken);
        }
    }

    auto clusterHits_after1 = getReconstructedHits(false);
    CHECK(clusterHits_after1.size() == 5);
    CHECK(clusterHits_after1[Detector_t::Type_t::CB] == 1835);
    CHECK(clusterHits_after1[Detector_t::Type_t::TAPS] == 359);
    CHECK(clusterHits_after1[Detector_t::Type_t::PID] == 321);
    CHECK(clusterHits_after1[Detector_t::Type_t::TAPSVeto] == 377);
    CHECK(clusterHits_after1[Detector_t::Type_t::EPT] == 3106);

    // ignore all elements in each detector type
    for(const auto& item : clusterHits_before) {
        auto det = ExpConfig::Setup::GetDetector(item.first);
        for(unsigned ch=0;ch<det->GetNChannels();ch++) {
            det->SetElementFlags(ch, Detector_t::ElementFlag_t::Broken);
        }
    }

    // nothing should be returned
    auto clusterHits_after2 = getReconstructedHits(false);
    CHECK(clusterHits_after2.empty());
}

void dotest_ignoredelements_raw_include() {
    auto clusterHits_before = getReconstructedHits(false);
    CHECK(clusterHits_before.size() == 5);
    CHECK(clusterHits_before[Detector_t::Type_t::CB] == 3759);
    CHECK(clusterHits_before[Detector_t::Type_t::TAPS] == 720);
    CHECK(clusterHits_before[Detector_t::Type_t::PID] == 633);
    CHECK(clusterHits_before[Detector_t::Type_t::TAPSVeto] == 956);
    CHECK(clusterHits_before[Detector_t::Type_t::EPT] == 6272);

    // ignore all elements in each detector type
    for(const auto& item : clusterHits_before) {
        auto det = ExpConfig::Setup::GetDetector(item.first);
        for(unsigned ch=0;ch<det->GetNChannels();ch++) {
            det->SetElementFlags(ch, Detector_t::ElementFlag_t::Broken);
        }
    }

    // nothing should have changed
    auto clusterHits_after2 = getReconstructedHits(false);
    CHECK(clusterHits_after2.size() == 5);
    CHECK(clusterHits_after2[Detector_t::Type_t::CB] == 3759);
    CHECK(clusterHits_after2[Detector_t::Type_t::TAPS] == 720);
    CHECK(clusterHits_after2[Detector_t::Type_t::PID] == 633);
    CHECK(clusterHits_after2[Detector_t::Type_t::TAPSVeto] == 956);
    CHECK(clusterHits_after2[Detector_t::Type_t::EPT] == 6272);
}

void dotest_ignoredelements_geant() {
    auto clusterHits_before = getReconstructedHits(true);
    CHECK(clusterHits_before.size() == 5);
    CHECK(clusterHits_before[Detector_t::Type_t::CB] == 1811);
    CHECK(clusterHits_before[Detector_t::Type_t::TAPS] == 675);
    CHECK(clusterHits_before[Detector_t::Type_t::PID] == 51);
    CHECK(clusterHits_before[Detector_t::Type_t::TAPSVeto] == 133);
    CHECK(clusterHits_before[Detector_t::Type_t::EPT] == 100);

    // ignore every second element in each detector type
    for(const auto& item : clusterHits_before) {
        auto det = ExpConfig::Setup::GetDetector(item.first);
        for(unsigned ch=0;ch<det->GetNChannels();ch++) {
            if(ch % 2 == 0)
                det->SetElementFlags(ch, Detector_t::ElementFlag_t::Broken);
        }
    }

    auto clusterHits_after1 = getReconstructedHits(true);
    CHECK(clusterHits_after1.size() == 5);
    CHECK(clusterHits_after1[Detector_t::Type_t::CB] == 896);
    CHECK(clusterHits_after1[Detector_t::Type_t::TAPS] == 329);
    CHECK(clusterHits_after1[Detector_t::Type_t::PID] == 25);
    CHECK(clusterHits_after1[Detector_t::Type_t::TAPSVeto] == 62);
    CHECK(clusterHits_after1[Detector_t::Type_t::EPT] == 54);

    // ignore all elements in each detector type
    for(const auto& item : clusterHits_before) {
        auto det = ExpConfig::Setup::GetDetector(item.first);
        for(unsigned ch=0;ch<det->GetNChannels();ch++) {
            det->SetElementFlags(ch, Detector_t::ElementFlag_t::Broken);
        }
    }

    // nothing should be returned
    auto clusterHits_after2 = getReconstructedHits(true);
    CHECK(clusterHits_after2.empty());
}

void dotest_ignoredelements_geant_include() {
    auto clusterHits_before = getReconstructedHits(true);
    CHECK(clusterHits_before.size() == 5);
    CHECK(clusterHits_before[Detector_t::Type_t::CB] == 1811);
    CHECK(clusterHits_before[Detector_t::Type_t::TAPS] == 675);
    CHECK(clusterHits_before[Detector_t::Type_t::PID] == 51);
    CHECK(clusterHits_before[Detector_t::Type_t::TAPSVeto] == 133);
    CHECK(clusterHits_before[Detector_t::Type_t::EPT] == 100);

    // ignore all elements in each detector type
    for(const auto& item : clusterHits_before) {
        auto det = ExpConfig::Setup::GetDetector(item.first);
        for(unsigned ch=0;ch<det->GetNChannels();ch++) {
            det->SetElementFlags(ch, Detector_t::ElementFlag_t::Broken);
        }
    }

    // nothing should have changed
    auto clusterHits_after2 = getReconstructedHits(true);
    CHECK(clusterHits_after2.size() == 5);
    CHECK(clusterHits_after2[Detector_t::Type_t::CB] == 1811);
    CHECK(clusterHits_after2[Detector_t::Type_t::TAPS] == 675);
    CHECK(clusterHits_after2[Detector_t::Type_t::PID] == 51);
    CHECK(clusterHits_after2[Detector_t::Type_t::TAPSVeto] == 133);
    CHECK(clusterHits_before[Detector_t::Type_t::EPT] == 100);
}