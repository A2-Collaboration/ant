#include "Reconstruct.h"

#include "HitMatching.h"
#include "Clustering.h"
#include "ApplyCalibrations.h"

#include "expconfig/ExpConfig.h"

#include "tree/TDataRecord.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"

#include <algorithm>
#include <iostream>

using namespace std;
using namespace ant;

Reconstruct::Reconstruct(const THeaderInfo &headerInfo)
{
  const auto& config = ExpConfig::Reconstruct::Get(headerInfo);
  calibrations = config->GetCalibrations();
  updateables = config->GetUpdateables();
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

  auto hint = sorted_clusterhits.begin();
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
        clusterhits.emplace_back(readhit->Channel, data);
      }
      else {
        // clusterhit with channel of readhit exists,
        // so append TClusterHitDatum's
        auto& data_ = it_clusterhit->Data;
        data_.insert(data_.end(), data.cbegin(), data.cend());
      }
    }

    // The trigger detector, for example, might only carry
    // reference hits, which are all not calibrated and thus
    // never fill anything in clusterhits
    if(clusterhits.empty())
      continue;

    hint = sorted_clusterhits.insert(hint, make_pair(detector, clusterhits));
  }


  // now we can do clustering on the sorted_clusterhits (by detector)

  return nullptr;
}
