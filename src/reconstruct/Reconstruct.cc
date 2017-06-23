#include "Reconstruct.h"

#include "Clustering.h"
#include "CandidateBuilder.h"
#include "UpdateableManager.h"

#include "expconfig/ExpConfig.h"

#include "tree/TEventData.h"

#include "base/std_ext/container.h"
#include "base/Logger.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::reconstruct;

Reconstruct::clustering_t Reconstruct::GetDefaultClustering()
{
    return std_ext::make_unique<Clustering_NextGen>();
}

Reconstruct::candidatebuilder_t Reconstruct::GetDefaultCandidateBuilder()
{
    /// \todo instead of using the full-blown CandidateBuilder here,
    /// it would be better to compose it out of smaller parts
    /// this might get important if using the MWPCs...
    try {
        return std_ext::make_unique<CandidateBuilder>();
    }
    catch(ExpConfig::ExceptionNoDetector e) {
        LOG(WARNING) << "Candidate builder could not be activated: " << e.what();
    }
    return nullptr;
}

template<typename List>
List getSortedHooks() {
    typename std::remove_const<List>::type hooks;
    using shared_ptr_t = typename List::value_type;
    using element_t = typename shared_ptr_t::element_type;
    for(const auto& hook : ExpConfig::Setup::Get().GetReconstructHooks()) {
        std_ext::AddToSharedPtrList<element_t>(hook, hooks);
    }
    return hooks;
}

Reconstruct::sorted_detectors_t Reconstruct::sorted_detectors_t::Build()
{
    sorted_detectors_t sorted_detectors;
    for(const auto& detector : ExpConfig::Setup::Get().GetDetectors()) {
        auto ret = sorted_detectors.insert(make_pair(detector->Type, detector));
        if(!ret.second) {
            throw Exception("Setup provided detector list with two detectors of same type");
        }
    }
    return sorted_detectors;
}

Reconstruct::Reconstruct(clustering_t clustering_, candidatebuilder_t candidatebuilder_) :
    includeIgnoredElements(ExpConfig::Setup::Get().GetIncludeIgnoredElements()),
    sorted_detectors(sorted_detectors_t::Build()),
    hooks_readhits(getSortedHooks<decltype(hooks_readhits)>()),
    hooks_clusterhits(getSortedHooks<decltype(hooks_clusterhits)>()),
    hooks_clusters(getSortedHooks<decltype(hooks_clusters)>()),
    hooks_eventdata(getSortedHooks<decltype(hooks_eventdata)>()),
    clustering(move(clustering_)),
    candidatebuilder(move(candidatebuilder_)),
    updateablemanager(std_ext::make_unique<UpdateableManager>(ExpConfig::Setup::Get().GetUpdateables()))
{
}

// implement the destructor here,
// makes forward declaration work properly
Reconstruct::~Reconstruct() = default;

void Reconstruct::DoReconstruct(TEventData& reconstructed) const
{
    // ignore empty events
    if(reconstructed.DetectorReadHits.empty())
        return;

    // update the updateables :)
    updateablemanager->UpdateParameters(reconstructed.ID);

    // apply the hooks for detector read hits (mostly calibrations),
    // note that this also changes the hits itself

    ApplyHooksToReadHits(reconstructed.DetectorReadHits);
    // the detectorReads are now calibrated as far as possible
    // one might return now and detectorRead is just calibrated...


    // do the hit matching, which builds the TClusterHit's
    // put into the AdaptorTClusterHit to track Energy/Timing information
    // for subsequent clustering
    sorted_bydetectortype_t<TClusterHit> sorted_clusterhits;
    BuildHits(sorted_clusterhits, reconstructed.TaggerHits);

    // apply hooks which modify clusterhits
    for(const auto& hook : hooks_clusterhits) {
        hook->ApplyTo(sorted_clusterhits);
    }

    // then build clusters (at least for calorimeters this is not trivial)
    sorted_clusters_t sorted_clusters;
    BuildClusters(move(sorted_clusterhits), sorted_clusters);

    // apply hooks which modify clusters
    for(const auto& hook : hooks_clusters) {
        hook->ApplyTo(sorted_clusters);
    }

    // do the candidate building (if available)
    if(candidatebuilder)
        candidatebuilder->Build(move(sorted_clusters),
                                reconstructed.Candidates, reconstructed.Clusters);

    // apply hooks which may modify the whole event
    for(const auto& hook : hooks_eventdata) {
        hook->ApplyTo(reconstructed);
    }

}

void Reconstruct::ApplyHooksToReadHits(std::vector<TDetectorReadHit>& detectorReadHits) const
{
    // categorize the hits by detector type
    // this is handy for all subsequent reconstruction steps
    // we need to use non-const references because calibrations
    // may change the content (use std::reference_wrapper to hold it in vector)
    sorted_readhits.clear();
    for(TDetectorReadHit& readhit : detectorReadHits) {
        sorted_readhits.add_item(readhit.DetectorType, readhit);
    }

    // apply calibration
    // this may change the given readhits
    for(const auto& hook : hooks_readhits) {
        hook->ApplyTo(sorted_readhits);
    }
}

void Reconstruct::BuildHits(sorted_bydetectortype_t<TClusterHit>& sorted_clusterhits,
        vector<TTaggerHit>& taggerhits) const
{
    auto insert_hint = sorted_clusterhits.cbegin();

    for(const auto& it_hit : sorted_readhits) {
        const Detector_t::Type_t detectortype = it_hit.first;
        const auto& readhits = it_hit.second;

        // find the detector instance for this type
        const auto& it_detector = sorted_detectors.find(detectortype);
        if(it_detector == sorted_detectors.end())
            continue;
        const detector_ptr_t& detector = it_detector->second;

        // for tagger detectors, we do not match the hits by channel at all
        if(detector.TaggerDetector != nullptr) {
            HandleTagger(detector.TaggerDetector, readhits, taggerhits);
            continue;
        }

        map<unsigned, TClusterHit> hits;

        for(const TDetectorReadHit& readhit : readhits) {
            if(!includeIgnoredElements && detector.Detector->IsIgnored(readhit.Channel))
                continue;

            // ignore uncalibrated items
            if(readhit.Values.empty())
                continue;


            auto& clusterhit = hits[readhit.Channel];
            // copy over all readhit info to clusterhit
            // For example, CB_TimeWalk needs all timings here!
            for(auto& v : readhit.Values)
                clusterhit.Data.emplace_back(readhit.ChannelType, v);
            clusterhit.Channel = readhit.Channel; // copy over the channel

            // set the energy or timing field (might stay NaN if not calibrated)
            // for multihit timing
            if(readhit.ChannelType == Channel_t::Type_t::Integral)
                clusterhit.Energy = readhit.Values.front().Calibrated;
            else if(readhit.ChannelType == Channel_t::Type_t::Timing)
                clusterhit.Time = readhit.Values.front().Calibrated;
        }

        TClusterHitList clusterhits;
        for(auto& it_hit : hits) {
            auto& hit = it_hit.second;

            // check for weird energies
            if(hit.IsSane() && hit.Energy<0) {
                // mostly PID/TAPS/TAPSVeto channels with there pedestal subtraction
                // cause negative energy entries, but that should be handled by
                // a meaningful raw threshold (if not disabled in setup)
                VLOG(7) << "Cluster Hit Energy " << hit.Energy << " MeV less than zero, ignoring. Det="
                        << Detector_t::ToString(detectortype) << " Ch=" << hit.Channel;
                hit.Energy = std_ext::NaN;
            }
            clusterhits.emplace_back(move(it_hit.second));
        }


        // The trigger or tagger detectors don't fill anything
        // so skip it
        if(clusterhits.empty())
            continue;

        // insert the clusterhits
        insert_hint =
                sorted_clusterhits.insert(insert_hint,
                                          make_pair(detectortype, move(clusterhits)));
    }
}

void Reconstruct::HandleTagger(const shared_ptr<TaggerDetector_t>& taggerdetector,
                               const std::vector<std::reference_wrapper<TDetectorReadHit> >& readhits,
                               std::vector<TTaggerHit>& taggerhits
                               ) const
{

    // gather electron hits by channel
    struct taggerhit_t {
        std::vector<TDetectorReadHit::Value_t> Timings;
        std::vector<TDetectorReadHit::Value_t> Energies;
    };
    map<unsigned, taggerhit_t > hits;

    for(const TDetectorReadHit& readhit : readhits) {
        if(!includeIgnoredElements && taggerdetector->IsIgnored(readhit.Channel))
            continue;

        // ignore uncalibrated items
        if(readhit.Values.empty())
            continue;

        auto& item = hits[readhit.Channel];
        if(readhit.ChannelType == Channel_t::Type_t::Timing) {
            std_ext::concatenate(item.Timings, readhit.Values);
        }
        else if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            std_ext::concatenate(item.Energies, readhit.Values);
        }
    }

    for(const auto& hit : hits) {
        const auto channel = hit.first;
        const auto& item = hit.second;
        // create a taggerhit from each timing for now
        /// \todo handle double hits here?
        /// \todo handle energies here better? (actually test with appropiate QDC run)
        const auto qdc_energy = item.Energies.empty() ? std_ext::NaN : item.Energies.front().Calibrated;
        for(const auto& timing : item.Timings) {
            taggerhits.emplace_back(channel,
                                    taggerdetector->GetPhotonEnergy(channel),
                                    timing.Calibrated,
                                    qdc_energy
                                    );
        }
    }
}

void Reconstruct::BuildClusters(
        const sorted_clusterhits_t& sorted_clusterhits,
        sorted_clusters_t& sorted_clusters) const
{
    auto insert_hint = sorted_clusters.begin();

    for(const auto& it_clusterhits : sorted_clusterhits) {
        const Detector_t::Type_t detectortype = it_clusterhits.first;
        const TClusterHitList& clusterhits = it_clusterhits.second;

        // find the detector instance for this type
        const auto& it_detector = sorted_detectors.find(detectortype);
        if(it_detector == sorted_detectors.end())
            continue;
        const detector_ptr_t& detector = it_detector->second;

        TClusterList clusters;

        // check if detector supports clustering
        if(detector.ClusterDetector != nullptr) {
            // yes, then hand over to clustering algorithm (if available)
            if(clustering)
                clustering->Build(*detector.ClusterDetector, clusterhits, clusters);
        }
        else {
            // in case of no clustering detector,
            // build simple "cluster" consisting of single TClusterHit
            for(const TClusterHit& hit : clusterhits) {

                // ignore hits with time and energy information
                if(!hit.IsSane())
                    continue;


                clusters.emplace_back(
                                          detector.Detector->GetPosition(hit.Channel),
                                          hit.Energy,
                                          hit.Time,
                                          detector.Detector->Type,
                                          hit.Channel,
                                          vector<TClusterHit>{hit}

                                      );

            }
        }

        // insert the clusters (if any)
        if(!clusters.empty()) {
            insert_hint =
                    sorted_clusters.insert(insert_hint,
                                           make_pair(detectortype, move(clusters)));
        }
    }
}



