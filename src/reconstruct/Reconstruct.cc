#include "Reconstruct.h"

#include "Clustering.h"

#include "expconfig/ExpConfig.h"

#include "tree/TDataRecord.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"

#include <algorithm>
#include <iostream>
#include <iterator>

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
    // ... they may be a cluster detector
    const auto& detector_cluster = dynamic_pointer_cast<ClusterDetector_t, Detector_t>(detector);
    if(detector_cluster != nullptr)
      detectors_cluster.push_back(detector_cluster);
  }

  /// \todo build the range list from updateables
}

unique_ptr<TEvent> Reconstruct::DoReconstruct(TDetectorRead& detectorRead)
{

  // categorize the hits by detector type
  // this is handy for all subsequent reconstruction steps
  map<Detector_t::Type_t, std::list< TDetectorReadHit* > > sorted_readhits;
  for(TDetectorReadHit& readhit : detectorRead.Hits) {
    const auto detector = static_cast<Detector_t::Type_t>(readhit.DetectorType);
    sorted_readhits[detector].push_back(addressof(readhit));
  }

  // apply calibration (this alters the given detectorRead!)
  for(const auto& calib : calibrations) {
    calib->ApplyTo(sorted_readhits);
  }

  // for debug purposes, dump out the detectorRead
  //cout << detectorRead << endl;

  // the detectorRead is now calibrated as far as possible
  // lets start the hit matching, which builds the TClusterHit's
  map<Detector_t::Type_t, std::list< TClusterHit > > sorted_clusterhits;

  auto sorted_clusterhits_hint = sorted_clusterhits.begin(); // use some hinting to speed up the map filling
  for(const auto& it_hit : sorted_readhits) {
    std::list<TClusterHit> clusterhits;
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
      const auto match_channel = [readhit] (const TClusterHit& hit) {
        return hit.Channel == readhit->Channel;
      };
      const auto it_clusterhit = find_if(clusterhits.begin(),
                                         clusterhits.end(),
                                         match_channel);
      if(it_clusterhit == clusterhits.end()) {
        // not found, create new TClusterHit from readhit
        clusterhits.emplace_back(readhit->Channel, move(data));
      }
      else {
        // clusterhit with channel of readhit already exists,
        // so append TClusterHitDatum's
        move(data.begin(), data.end(),
             back_inserter(it_clusterhit->Data));
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


  // now we can do clustering on the sorted_clusterhits (sorted by detector)
  // we build some wrapper structure to hold necessary information for clustering





  return nullptr;
}
