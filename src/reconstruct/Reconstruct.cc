#include "Reconstruct.h"

#include "Clustering.h"
#include "CandidateBuilder.h"
#include "UpdateableManager.h"

#include "expconfig/ExpConfig.h"

#include "tree/TDataRecord.h"
#include "tree/TDetectorRead.h"
#include "tree/THeaderInfo.h"
#include "tree/TEvent.h"

#include "base/Logger.h"
#include "base/std_ext.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>

using namespace std;
using namespace ant;
using namespace ant::reconstruct;

Reconstruct::Reconstruct(const THeaderInfo& headerInfo)
{
    const auto& config = ExpConfig::Reconstruct::Get(headerInfo);

    // calibrations and detector may change their parameters
    // during reconstruct, so we scan for those items
    shared_ptr_list<Updateable_traits> updateables;

    // hooks are usually calibrations, which may also be updateable
    const shared_ptr_list<ReconstructHook::Base>& hooks = config->GetReconstructHooks();
    for(const auto& hook : hooks) {
        std_ext::AddToSharedPtrList<ReconstructHook::DetectorReadHits, ReconstructHook::Base>
                (hook, hooks_readhits);
        std_ext::AddToSharedPtrList<ReconstructHook::Clusters, ReconstructHook::Base>
                (hook, hooks_clusters);
        std_ext::AddToSharedPtrList<Updateable_traits, ReconstructHook::Base>
                (hook, updateables);
    }

    // detectors serve different purposes, ...
    const shared_ptr_list<Detector_t>& detectors = config->GetDetectors();
    for(const auto& detector : detectors) {
        // ... they may be updateable (think of the PID Phi angles)
        std_ext::AddToSharedPtrList<Updateable_traits, Detector_t>
                (detector, updateables);
        // ... but also are needed in DoReconstruct
        auto ret = sorted_detectors.insert(make_pair(detector->Type, detector));
        if(!ret.second) {
            throw runtime_error("Reconstruct config provided detector list with two detectors of same type");
        }
    }

    // init clustering
    clustering = std_ext::make_unique<Clustering>(config);

    // init the candidate builder
    /// \todo Make use of different candidate builders maybe?
    candidatebuilder = std_ext::make_unique<CandidateBuilder>(sorted_detectors);

    // init the updateable manager
    updateablemanager = std_ext::make_unique<UpdateableManager>(
                            headerInfo.ID, move(updateables)
                            );
}

// implement the destructor here,
// makes forward declaration work properly
Reconstruct::~Reconstruct() {}

unique_ptr<TEvent> Reconstruct::DoReconstruct(TDetectorRead& detectorRead)
{
    // update the updateables :)
    updateablemanager->UpdateParameters(detectorRead.ID);

    // apply the hooks for detector read hits (mostly calibrations),
    // note that this also changes the detectorRead

    ApplyHooksToReadHits(detectorRead);
    // the detectorRead is now calibrated as far as possible
    // one might return now and detectorRead is just calibrated...

    // for debug purposes, dump out the detectorRead
    //cout << detectorRead << endl;

    // already create the event here, since Tagger
    // doesn't need hit matching and thus can be filled already
    // in BuildHits (see below)
    auto event = std_ext::make_unique<TEvent>(detectorRead.ID);

    // do the hit matching, which builds the TClusterHit's
    // we also extract the energy, which is always defined as a
    // single value with type Channel_t::Type_t
    sorted_bydetectortype_t<AdaptorTClusterHit> sorted_clusterhits;
    BuildHits(sorted_clusterhits, event->Tagger);

    // then build clusters (at least for calorimeters this is not trivial)
    sorted_bydetectortype_t<TCluster> sorted_clusters;
    BuildClusters(move(sorted_clusterhits), sorted_clusters, event->InsaneClusters);

    // apply hooks which modify clusters
    for(const auto& hook : hooks_clusters) {
        hook->ApplyTo(sorted_clusters);
    }

    // finally, do the candidate building
    candidatebuilder->Build(move(sorted_clusters), event->Candidates);

    // uncomment for debug purposes
    //cout << *event << endl;

    return event;
}

void Reconstruct::ApplyHooksToReadHits(
        TDetectorRead& detectorRead)
{
    // categorize the hits by detector type
    // this is handy for all subsequent reconstruction steps
    // we need to use non-const pointers because calibrations
    // may change the content
    sorted_readhits.clear();
    for(TDetectorReadHit& readhit : detectorRead.Hits) {
        sorted_readhits.add_item(readhit.GetDetectorType(), addressof(readhit));
    }

    // apply calibration
    // this may change the given detectorRead,
    // and even provide use with extrahits to be added to it
    list<TDetectorReadHit> extrahits;
    for(const auto& hook : hooks_readhits) {
        hook->ApplyTo(sorted_readhits, extrahits);
    }

    // adding extrahits to the detectorRead invalidates the sorted_readhits pointer map
    /// \todo make this more efficient
    if(!extrahits.empty()) {
        detectorRead.Hits.insert(detectorRead.Hits.end(),
                                 make_move_iterator(extrahits.begin()),
                                 make_move_iterator(extrahits.end()));
        sorted_readhits.clear();
        for(TDetectorReadHit& readhit : detectorRead.Hits) {
            sorted_readhits.add_item(readhit.GetDetectorType(), addressof(readhit));
        }
    }
}

void Reconstruct::BuildHits(sorted_bydetectortype_t<AdaptorTClusterHit>& sorted_clusterhits,
        TTagger& event_tagger)
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
        const shared_ptr<Detector_t>& detector = it_detector->second;

        // the tagging devices are excluded from further hit matching and clustering
        const shared_ptr<TaggerDetector_t>& taggerdetector
                = dynamic_pointer_cast<TaggerDetector_t>(detector);

        // for tagger detectors, we do not match the hits by channel at all
        if(taggerdetector != nullptr) {
            HandleTagger(taggerdetector, readhits, event_tagger);
            continue;
        }

        for(const TDetectorReadHit* readhit : readhits) {
            // ignore uncalibrated items
            if(readhit->Values.empty())
                continue;



            // transform the data from readhit into TClusterHitDatum's
            vector<TClusterHitDatum> data(readhit->Values.size());
            auto do_transform = [readhit] (double value) {
                return TClusterHitDatum(readhit->GetChannelType(), value);
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
                               const std::vector<TDetectorReadHit*>& readhits,
                               TTagger& event_tagger
                               )
{
    /// \todo add Moeller/PairSpec information here

    for(const TDetectorReadHit* readhit : readhits) {
        // ignore uncalibrated items
        if(readhit->Values.empty())
            continue;

        switch(readhit->GetChannelType()) {

        case Channel_t::Type_t::Timing: {
            // but we add each TClusterHitDatum as some single
            // TClusterHit representing an electron with a timing
            /// \todo implement tagger double-hit decoding?
            for(const double timing : readhit->Values) {
                event_tagger.Hits.emplace_back(
                            taggerdetector->GetPhotonEnergy(readhit->Channel),
                            TKeyValue<double>(readhit->Channel, timing)
                            );
            }
            break;
        }

        case Channel_t::Type_t::Scaler: {
            // scalers are already calibrated to frequencies
            for(const double frequency : readhit->Values) {
                event_tagger.Scalers.emplace_back(readhit->Channel, frequency);
            }
            break;
        }

        case Channel_t::Type_t::Integral: {
            // tagger energies are rarely present, only for special fastbus QDC runs
            /// \todo actually test this part of the code
            for(const double energy : readhit->Values) {
                event_tagger.Energies.emplace_back(readhit->Channel, energy);
            }
            break;
        }

        default:
            break;
        }

    }
}

void Reconstruct::BuildClusters(
        sorted_bydetectortype_t<AdaptorTClusterHit>&& sorted_clusterhits,
        sorted_bydetectortype_t<TCluster>& sorted_clusters,
        std::vector<TCluster>& insane_clusters)
{
    auto insert_hint = sorted_clusters.begin();

    for(const auto& it_clusterhits : sorted_clusterhits) {
        const Detector_t::Type_t detectortype = it_clusterhits.first;
        const list<AdaptorTClusterHit>& clusterhits = it_clusterhits.second;

        // find the detector instance for this type
        const auto& it_detector = sorted_detectors.find(detectortype);
        if(it_detector == sorted_detectors.end())
            continue;
        const shared_ptr<Detector_t>& detector = it_detector->second;

        list<TCluster> clusters;

        // check if detector can do clustering,
        const shared_ptr<ClusterDetector_t>& clusterdetector
                = dynamic_pointer_cast<ClusterDetector_t>(detector);

        if(clusterdetector != nullptr) {
            // yes, then hand over to clustering algorithm
            clustering->Build(clusterdetector, clusterhits, clusters);
        }
        else {
            // in case of no clustering detector,
            // build simple "cluster" consisting of single TClusterHit
            for(const AdaptorTClusterHit& clusterhit : clusterhits) {
                const auto& hit = clusterhit.Hit;

                clusters.emplace_back(
                            detector->GetPosition(hit->Channel),
                            clusterhit.Energy,
                            clusterhit.Time,
                            detector->Type,
                            hit->Channel,
                            vector<TClusterHit>{*hit}
                            );

            }
        }

        auto c = clusters.begin();
        while (c!=clusters.end()) {
            if(! c->isSane()) {
                insane_clusters.emplace_back(std::move(*c));
                c = clusters.erase(c);
            } else {
                ++c;
            }
        }

        if(clusters.empty())
            continue;

        // insert the clusters
        insert_hint =
                sorted_clusters.insert(insert_hint,
                                       make_pair(detectortype, move(clusters)));
    }
}

