#include "Reconstruct.h"

#include "Clustering.h"
#include "CandidateBuilder.h"

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

    // calibrations are the natural Updateable_traits objects,
    // but they could be constant (in case of a very simple calibration)
    calibrations = config->GetCalibrations();
    for(const auto& calibration : calibrations) {
        const auto& ptr = dynamic_pointer_cast<Updateable_traits, CalibrationApply_traits>(calibration);
        if(ptr != nullptr)
            updateables.push_back(ptr);
    }

    // detectors serve different purposes, ...
    auto detectors = config->GetDetectors();
    for(const auto& detector : detectors) {
        // ... they may be updateable (think of the PID Phi angles)
        const auto& updateable = dynamic_pointer_cast<Updateable_traits, Detector_t>(detector);
        if(updateable != nullptr)
            updateables.push_back(updateable);
        // ... but also are needed in DoReconstruct
        auto ret = sorted_detectors.insert(make_pair(detector->Type, detector));
        if(!ret.second) {
            throw runtime_error("Reconstruct config provided detector list with two detectors of same type");
        }
    }

    // init clustering
    clustering = std_ext::make_unique<Clustering>(config);

    // init the trackbuilder
    /// \todo Make use of different TrackBuilders maybe?
    candidatebuilder = std_ext::make_unique<CandidateBuilder>(sorted_detectors);

    // ask each updateable for its update points and
    // build the list of changePoints

    map<TID, shared_ptr_list<Updateable_traits> > sorted_updateables;
    for(const shared_ptr<Updateable_traits>& updateable : updateables) {
        for(const TID& changePoint : updateable->GetChangePoints()) {
            // ignore too early change points
            /// \todo maybe headerInfo timepoint lies in between change points?
            if(changePoint < headerInfo.ID)
                continue;
            sorted_updateables[changePoint].emplace_back(move(updateable));
        }
    }

    // we rely on the fact that the map is sorted,
    // and simply convert it to a list
    for(auto it_updateable : sorted_updateables) {
        changePoints.emplace_back(it_updateable);
    }
}

unique_ptr<TEvent> Reconstruct::DoReconstruct(TDetectorRead& detectorRead)
{
    // update the updateables :)
    if(!changePoints.empty() && changePoints.front().first <= detectorRead.ID) {
        unsigned nUpdateables = 0;
        for(const auto& updateable : changePoints.front().second) {
            updateable->Update(changePoints.front().first);
            nUpdateables++;
        }

        VLOG(7) << "Updated parameters for " << nUpdateables
                << " calibrations and detectors at ID=" << detectorRead.ID;
        changePoints.pop_front(); // go to next change point
    }

    // apply the calibrations,
    // note that this also changes the detectorRead
    CalibrationApply_traits::readhits_t sorted_readhits;
    ApplyCalibrations(detectorRead, sorted_readhits);

    // for debug purposes, dump out the detectorRead
    //cout << detectorRead << endl;

    // already create the event here, since Tagger
    // doesn't need hit matching and thus can be filled already
    // in BuildHits (see below)
    auto event = std_ext::make_unique<TEvent>(detectorRead.ID);

    // the detectorRead is now calibrated as far as possible
    // lets start the hit matching, which builds the TClusterHit's
    // we also extract the energy, which is always defined as a
    // single value with type Channel_t::Type_t
    sorted_bydetectortype_t<AdaptorTClusterHit> sorted_clusterhits;
    BuildHits(move(sorted_readhits), sorted_clusterhits, event->Tagger);

    // then build clusters (at least for calorimeters this is not trivial)
    sorted_bydetectortype_t<TCluster> sorted_clusters;
    BuildClusters(move(sorted_clusterhits), sorted_clusters, event->InsaneClusters);


    // finally, do the candidate building
    candidatebuilder->Build(move(sorted_clusters), event->Candidates);

    // uncomment for debug purposes
    //cout << *event << endl;

    return event;
}

Reconstruct::~Reconstruct()
{

}

void Reconstruct::ApplyCalibrations(TDetectorRead& detectorRead,
                                    sorted_bydetectortype_t<TDetectorReadHit*>& sorted_readhits)
{
    // categorize the hits by detector type
    // this is handy for all subsequent reconstruction steps
    // we need to use non-const pointers because calibrations
    // may change the content
    for(TDetectorReadHit& readhit : detectorRead.Hits) {
        sorted_readhits[readhit.GetDetectorType()].push_back(addressof(readhit));
    }

    // apply calibration
    // this may change the given detectorRead,
    // and even provide use with extrahits to be added to it
    list<TDetectorReadHit> extrahits;
    for(const auto& calib : calibrations) {
        calib->ApplyTo(sorted_readhits, extrahits);
    }

    // adding extrahits to the detectorRead invalidates the sorted_readhits pointer map
    /// \todo make this more efficient
    if(!extrahits.empty()) {
        detectorRead.Hits.insert(detectorRead.Hits.end(),
                                 make_move_iterator(extrahits.begin()),
                                 make_move_iterator(extrahits.end()));
        sorted_readhits.clear();
        for(TDetectorReadHit& readhit : detectorRead.Hits) {
            sorted_readhits[readhit.GetDetectorType()].push_back(addressof(readhit));
        }
    }
}

void Reconstruct::BuildHits(
        sorted_bydetectortype_t<TDetectorReadHit *>&& sorted_readhits,
        sorted_bydetectortype_t<AdaptorTClusterHit>& sorted_clusterhits,
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

        for(const TDetectorReadHit* readhit : readhits) {
            // ignore uncalibrated items
            if(readhit->Values.empty())
                continue;

            // for tagger detectors, we do not match the hits by channel at all
            if(taggerdetector != nullptr) {
                /// \todo handle Integral and Scaler type information here
                if(readhit->GetChannelType() != Channel_t::Type_t::Timing)
                    continue;
                // but we add each TClusterHitDatum as some single
                // TClusterHit representing an electron with a timing
                /// \todo implement tagger double-hit decoding?
                for(const double timing : readhit->Values) {
                    event_tagger.Hits.emplace_back(
                                taggerdetector->GetPhotonEnergy(readhit->Channel),
                                TKeyValue<double>(readhit->Channel, timing)
                                );
                }

                /// \todo add Moeller/PairSpec information here
                continue;
            }

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

        // insert the clusters
        insert_hint =
                sorted_clusters.insert(insert_hint,
                                       make_pair(detectortype, move(clusters)));
    }
}

