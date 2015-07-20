#include "Reconstruct.h"

#include "Clustering.h"

#include "expconfig/ExpConfig.h"

#include "tree/TDataRecord.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"

#include "base/std_ext.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>

using namespace std;
using namespace ant;
using namespace ant::reconstruct;

Reconstruct::Reconstruct(const THeaderInfo &headerInfo)
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
    sorted_detectors[detector->Type] = detector;
  }

  /// \todo build the range list from updateables
}

unique_ptr<TEvent> Reconstruct::DoReconstruct(TDetectorRead& detectorRead)
{

  // categorize the hits by detector type
  // this is handy for all subsequent reconstruction steps
  map<Detector_t::Type_t, list< TDetectorReadHit* > > sorted_readhits;
  for(TDetectorReadHit& readhit : detectorRead.Hits) {
    const auto detector = static_cast<Detector_t::Type_t>(readhit.DetectorType);
    sorted_readhits[detector].push_back(addressof(readhit));
  }

  // apply calibration (this may change the given detectorRead!)
  for(const auto& calib : calibrations) {
    calib->ApplyTo(sorted_readhits);
  }

  // for debug purposes, dump out the detectorRead
  //cout << detectorRead << endl;

  // the detectorRead is now calibrated as far as possible
  // lets start the hit matching, which builds the TClusterHit's
  // we also extract the energy, which is always defined as a
  // single value with type Channel_t::Type_t

  struct HitWithEnergy_t {
    TClusterHit Hit;
    double Energy = numeric_limits<double>::quiet_NaN();
    void MaybeSetEnergy(const TDetectorReadHit* readhit) {
      if(readhit->GetChannelType() != Channel_t::Type_t::Integral)
        return;
      if(readhit->Values.size() != 1)
        return;
      Energy = readhit->Values[0];
    }
    HitWithEnergy_t(const TDetectorReadHit* readhit,
                    const vector<TClusterHitDatum>&& data) :
      Hit(readhit->Channel, data)
    {
      MaybeSetEnergy(readhit);
    }
  };

  map<Detector_t::Type_t, list< HitWithEnergy_t > > sorted_clusterhits;
  auto sorted_clusterhits_hint = sorted_clusterhits.cbegin();

  for(const auto& it_hit : sorted_readhits) {
    list<HitWithEnergy_t> clusterhits;
    const Detector_t::Type_t detector = it_hit.first;
    const auto& readhits = it_hit.second;

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

      // search for a TClusterHit with same channel
      const auto match_channel = [readhit] (const HitWithEnergy_t& hit) {
        return hit.Hit.Channel == readhit->Channel;
      };
      const auto it_clusterhit = find_if(clusterhits.begin(),
                                         clusterhits.end(),
                                         match_channel);
      if(it_clusterhit == clusterhits.end()) {
        // not found, create new TClusterHit from readhit
        clusterhits.emplace_back(readhit, move(data));
      }
      else {
        // clusterhit with channel of readhit already exists,
        // so append TClusterHitDatum's and set energy
        it_clusterhit->MaybeSetEnergy(readhit);
        move(data.begin(), data.end(),
             back_inserter(it_clusterhit->Hit.Data));
      }
    }

    // The trigger detector, for example, might only carry
    // reference hits, which are all not calibrated and thus
    // never fill anything in clusterhits
    if(clusterhits.empty())
      continue;

    // insert the clusterhits
    sorted_clusterhits_hint =
        sorted_clusterhits.insert(sorted_clusterhits_hint,
                                  make_pair(detector, move(clusterhits)));
  }

  // already create the event here, since TaggerHits
  // don't need tracking and thus can be filled already
  auto event = std_ext::make_unique<TEvent>(detectorRead.ID);

  map<Detector_t::Type_t, list< TCluster > > sorted_clusters;
  auto sorted_clusters_hint = sorted_clusters.begin();

  for(const auto& it_clusterhits : sorted_clusterhits) {
    const Detector_t::Type_t detectortype = it_clusterhits.first;
    const list<HitWithEnergy_t>& clusterhits = it_clusterhits.second;

    // do we have a detector instance for this type?
    const auto& it_detector = sorted_detectors.find(detectortype);
    if(it_detector == sorted_detectors.end())
      continue;
    const shared_ptr<Detector_t>& detector = it_detector->second;

    // we handle stuff a bit differently for
    // each detector from now on: tagger, cluster, everything else
    list<TCluster> clusters;

    // first the tagging device, which is excluded from track matching
    const shared_ptr<TaggerDetector_t>& taggerdetector
        = dynamic_pointer_cast<TaggerDetector_t>(detector);
    if(taggerdetector != nullptr) {
      // one might do some double-hit decoding here...?
      /// \todo handle the TaggerHit stuff here, maybe include PairSpec and Moeller?
      continue;
    }

    // check if detector can do clustering,
    const shared_ptr<ClusterDetector_t>& clusterdetector
        = dynamic_pointer_cast<ClusterDetector_t>(detector);

    if(clusterdetector != nullptr) {
      // clustering detector, so we need additional information
      // to build the crystals_t
      list<clustering::crystal_t> crystals;
      for(const HitWithEnergy_t& clusterhit : clusterhits) {
        const TClusterHit& hit = clusterhit.Hit;
        if(!isfinite(clusterhit.Energy))
          continue;
        crystals.emplace_back(
              clusterhit.Energy,
              clusterdetector->GetClusterElement(hit.Channel),
              addressof(hit)
              );

      }
      // do the clustering
      vector< vector< clustering::crystal_t> > crystal_clusters;
      do_clustering(crystals, crystal_clusters);

      // now calculate some cluster properties,
      // and create TCluster out of it
      for(const auto& cluster : crystal_clusters) {
        const double cluster_energy = clustering::calc_total_energy(cluster);
        /// \todo already cut low-energy clusters here?
        TVector3 weightedPosition(0,0,0);
        double weightedSum = 0;
        std::vector<TClusterHit> clusterhits;
        clusterhits.reserve(cluster.size());
        for(const clustering::crystal_t& crystal : cluster) {
          double wgtE = clustering::calc_energy_weight(crystal.Energy, cluster_energy);
          weightedPosition += crystal.Element->Position * wgtE;
          weightedSum += wgtE;
          clusterhits.emplace_back(*crystal.Hit);
        }
        weightedPosition *= 1.0/weightedSum;
        clusters.emplace_back(
              weightedPosition,
              cluster_energy,
              detector->Type,
              clusterhits
              );
      }
    }
    else {
      // in case of no clustering detector,
      // build "cluster" consisting of single TClusterHit
      for(const HitWithEnergy_t& clusterhit : clusterhits) {
        const TClusterHit& hit = clusterhit.Hit;
        clusters.emplace_back(
              detector->GetPosition(hit.Channel),
              clusterhit.Energy,
              detector->Type,
              vector<TClusterHit>{hit}
              );
      }
    }

    // insert the clusters
    sorted_clusters_hint =
        sorted_clusters.insert(sorted_clusters_hint,
                               make_pair(detectortype, move(clusters)));
  }

  // super-stupid track matching for now
  // each track contains just one cluster

  for(const auto& it_clusters : sorted_clusters) {
    //const Detector_t::Type_t detectortype = it_clusters.first;
    const list<TCluster>& clusters = it_clusters.second;

    for(const TCluster& cluster : clusters) {
      event->Tracks.emplace_back(
            cluster.Energy,
            0, // time unknown...
            cluster.Position.Theta(),
            cluster.Position.Phi(),
            std::vector<TCluster>{cluster}
            );
    }
  }

  // uncomment for debug purposes
  //cout << *event << endl;

  return event;
}
