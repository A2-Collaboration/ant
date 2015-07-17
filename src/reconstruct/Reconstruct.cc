#include "Reconstruct.h"

#include "Clustering.h"

#include "expconfig/ExpConfig.h"

#include "tree/TDataRecord.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>

using namespace std;
using namespace ant;

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
    void SetEnergy(const TDetectorReadHit* readhit) {
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
      SetEnergy(readhit);
    }
  };

  map<Detector_t::Type_t, list< HitWithEnergy_t > > sorted_clusterhits;
  auto sorted_clusterhits_hint = sorted_clusterhits.begin();

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
        if(readhit->GetChannelType() == Channel_t::Type_t::Integral)
          clusterhits.emplace_back(readhit, move(data));
      }
      else {
        // clusterhit with channel of readhit already exists,
        // so append TClusterHitDatum's
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
  auto event = make_shared<TEvent>(detectorRead.ID);

  map<Detector_t::Type_t, list< TCluster > > sorted_clusters;
  auto sorted_clusters_hint = sorted_clusters.begin();

  for(const auto& it_clusterhits : sorted_clusterhits) {
    list<TCluster> clusters;
    const Detector_t::Type_t detectortype = it_clusterhits.first;
    const list<HitWithEnergy_t>& clusterhits = it_clusterhits.second;

    // do we have a detector instance for this type?
    const auto& it_detector = sorted_detectors.find(detectortype);
    if(it_detector == sorted_detectors.end())
      continue;
    const shared_ptr<Detector_t>& detector = it_detector->second;

    // we handle stuff a bit differently for
    // each detector from now on: tagger, cluster, everything else

    // first the tagging device, which is excluded from track matching
    if(detectortype == Detector_t::Type_t::Tagger) {
      // one might do some double-hit decoding here...?
      /// \todo handle the TaggerHit stuff here, maybe include PairSpec and Moeller?
      continue;
    }

    // check if detector can do clustering,
    const shared_ptr<ClusterDetector_t>& clusterdetector
        = dynamic_pointer_cast<ClusterDetector_t>(detector);

    if(clusterdetector != nullptr) {
      // clustering detector, we need more than energy and position

    }
    else {
      // in case of no clustering,
      // build cluster for each hit
      for(const HitWithEnergy_t& clusterhit : clusterhits) {

        clusters.emplace_back(TCluster());
      }
    }




    // insert the clusterhits
    sorted_clusters_hint =
        sorted_clusters.insert(sorted_clusters_hint,
                               make_pair(detectortype, move(clusters)));
  }



  return nullptr;
}
