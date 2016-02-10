#include "Reconstruct.h"

#include "Clustering.h"
#include "CandidateBuilder.h"
#include "UpdateableManager.h"

#include "expconfig/ExpConfig.h"

#include "tree/TEventData.h"

#include "base/std_ext/vector.h"
#include "base/Logger.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::reconstruct;

Reconstruct::Reconstruct() {}

// implement the destructor here,
// makes forward declaration work properly
Reconstruct::~Reconstruct() {}

void Reconstruct::Initialize(const TID& tid)
{
    const auto& config = ExpConfig::Reconstruct::Get(tid);

    // hooks are usually calibrations, which may also be updateable
    const shared_ptr_list<ReconstructHook::Base>& hooks = config->GetReconstructHooks();
    for(const auto& hook : hooks) {
        std_ext::AddToSharedPtrList<ReconstructHook::DetectorReadHits, ReconstructHook::Base>
                (hook, hooks_readhits);
        std_ext::AddToSharedPtrList<ReconstructHook::ClusterHits, ReconstructHook::Base>
                (hook, hooks_clusterhits);
        std_ext::AddToSharedPtrList<ReconstructHook::Clusters, ReconstructHook::Base>
                (hook, hooks_clusters);
        std_ext::AddToSharedPtrList<ReconstructHook::EventData, ReconstructHook::Base>
                (hook, hooks_eventdata);
    }

    // put the detectors in a map for convenient access
    const shared_ptr_list<Detector_t>& detectors = config->GetDetectors();
    for(const auto& detector : detectors) {
        auto ret = sorted_detectors.insert(make_pair(detector->Type, detector));
        if(!ret.second) {
            throw Exception("Reconstruct config provided detector list with two detectors of same type");
        }
    }

    // init clustering
    clustering = std_ext::make_unique<Clustering>(config);

    // init the candidate builder
    /// \todo Make use of different candidate builders maybe?
    candidatebuilder = std_ext::make_unique<CandidateBuilder>(sorted_detectors, config);

    // init the updateable manager
    updateablemanager = std_ext::make_unique<UpdateableManager>(
                            tid, config->GetUpdateables()
                            );

    initialized = true;
}

void Reconstruct::DoReconstruct(TEventData& reconstructed)
{
    // ignore empty events
    if(reconstructed.DetectorReadHits.empty())
        return;

    if(!initialized) {
        Initialize(reconstructed.ID);
    }

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
    sorted_bydetectortype_t<AdaptorTClusterHit> sorted_clusterhits;
    BuildHits(sorted_clusterhits, reconstructed.TaggerHits);

    // apply hooks which modify clusterhits
    for(const auto& hook : hooks_clusterhits) {
        hook->ApplyTo(sorted_clusterhits);
    }

    // then build clusters (at least for calorimeters this is not trivial)
    sorted_bydetectortype_t<TClusterPtr> sorted_clusters;
    BuildClusters(move(sorted_clusterhits), sorted_clusters);

    // apply hooks which modify clusters
    for(const auto& hook : hooks_clusters) {
        hook->ApplyTo(sorted_clusters);
    }

    // do the candidate building
    candidatebuilder->Build(move(sorted_clusters),
                            reconstructed.Candidates, reconstructed.Clusters);

    // apply hooks which may modify the whole event
    for(const auto& hook : hooks_eventdata) {
        hook->ApplyTo(reconstructed);
    }

}

void Reconstruct::ApplyHooksToReadHits(
        std::vector<TDetectorReadHit>& detectorReadHits)
{
    // categorize the hits by detector type
    // this is handy for all subsequent reconstruction steps
    // we need to use non-const pointers because calibrations
    // may change the content
    sorted_readhits.clear();
    for(TDetectorReadHit& readhit : detectorReadHits) {
        sorted_readhits.add_item(readhit.DetectorType, addressof(readhit));
    }

    // apply calibration
    // this may change the given readhits
    for(const auto& hook : hooks_readhits) {
        hook->ApplyTo(sorted_readhits);
    }
}

void Reconstruct::BuildHits(sorted_bydetectortype_t<AdaptorTClusterHit>& sorted_clusterhits,
        vector<TTaggerHit>& taggerhits)
{
    auto insert_hint = sorted_clusterhits.cbegin();

    for(const auto& it_hit : sorted_readhits) {
        list<AdaptorTClusterHit> clusterhits;
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

        for(const TDetectorReadHit* readhit : readhits) {
            // ignore uncalibrated items
            if(readhit->Values.empty())
                continue;



            // transform the data from readhit into TClusterHitDatum's
            vector<TClusterHitDatum> data(readhit->Values.size());
            const Channel_t::Type_t channeltype = readhit->ChannelType;
            auto do_transform = [channeltype] (double value) {
                return TClusterHitDatum(channeltype, value);
            };
            transform(readhit->Values.cbegin(), readhit->Values.cend(),
                      data.begin(), do_transform);

            // for non-tagger detectors, search for a TClusterHit with same channel
            const auto match_channel = [readhit] (const AdaptorTClusterHit& hit) {
                return hit.Hit->Channel == readhit->Channel;
            };
            const auto it_clusterhit = find_if(clusterhits.begin(),
                                               clusterhits.end(),
                                               match_channel);
            if(it_clusterhit == clusterhits.end()) {
                // not found, create new HitWithEnergy from readhit
                clusterhits.emplace_back(readhit, move(data));
            }
            else {
                // clusterhit with channel of readhit already exists,
                // so append TClusterHitDatum's and try to set properties
                it_clusterhit->SetFields(readhit);
                move(data.begin(), data.end(),
                     back_inserter(it_clusterhit->Hit->Data));
            }
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
                               const vector<TDetectorReadHit*>& readhits,
                               std::vector<TTaggerHit>& taggerhits
                               )
{

    // gather electron hits by channel
    struct taggerhit_t {
        std::vector<double> Timings;
        std::vector<double> Energies;
    };
    map<unsigned, taggerhit_t > hits;

    for(const TDetectorReadHit* readhit : readhits) {
        // ignore uncalibrated items
        if(readhit->Values.empty())
            continue;

        auto& item = hits[readhit->Channel];
        if(readhit->ChannelType == Channel_t::Type_t::Timing) {
            std_ext::concatenate(item.Timings, readhit->Values);
        }
        else if(readhit->ChannelType == Channel_t::Type_t::Integral) {
            std_ext::concatenate(item.Energies, readhit->Values);
        }
    }

    for(const auto& hit : hits) {
        const auto channel = hit.first;
        const auto& item = hit.second;
        // create a taggerhit from each timing for now
        /// \todo handle double hits here?
        /// \todo handle energies here better? (actually test with appropiate QDC run)
        const auto qdc_energy = item.Energies.empty() ? std_ext::NaN : item.Energies.front();
        for(const auto timing : item.Timings) {
            taggerhits.emplace_back(channel,
                                    taggerdetector->GetPhotonEnergy(channel),
                                    timing,
                                    qdc_energy
                                    );
        }
    }
}

void Reconstruct::BuildClusters(sorted_bydetectortype_t<AdaptorTClusterHit>&& sorted_clusterhits,
        sorted_bydetectortype_t<TClusterPtr>& sorted_clusters)
{
    auto insert_hint = sorted_clusters.begin();

    for(const auto& it_clusterhits : sorted_clusterhits) {
        const Detector_t::Type_t detectortype = it_clusterhits.first;
        const list<AdaptorTClusterHit>& clusterhits = it_clusterhits.second;

        // find the detector instance for this type
        const auto& it_detector = sorted_detectors.find(detectortype);
        if(it_detector == sorted_detectors.end())
            continue;
        const detector_ptr_t& detector = it_detector->second;

        list<TClusterPtr> clusters;

        // check if detector supports clustering
        if(detector.ClusterDetector != nullptr) {
            // yes, then hand over to clustering algorithm
            clustering->Build(detector.ClusterDetector, clusterhits, clusters);
        }
        else {
            // in case of no clustering detector,
            // build simple "cluster" consisting of single TClusterHit
            for(const AdaptorTClusterHit& clusterhit : clusterhits) {

                // ignore hits with time and energy information
                if(!isfinite(clusterhit.Energy) || !isfinite(clusterhit.Time))
                    continue;

                const auto& hit = clusterhit.Hit;

                clusters.emplace_back(make_shared<TCluster>(
                                          detector.Detector->GetPosition(hit->Channel),
                                          clusterhit.Energy,
                                          clusterhit.Time,
                                          detector.Detector->Type,
                                          hit->Channel,
                                          vector<TClusterHit>{*hit}
                                          )
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

