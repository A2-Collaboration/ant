#include "Reconstruct.h"

#include "HitMatching.h"
#include "Clustering.h"
#include "ApplyCalibrations.h"

#include "expconfig/ExpConfig.h"

#include "tree/TDataRecord.h"
#include "tree/TDetectorRead.h"

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
  map<Detector_t::Type_t, std::list< TDetectorReadHit* > > sorted_hits;
  for(TDetectorReadHit& hit : detectorRead.Hits) {
    const auto detector = static_cast<Detector_t::Type_t>(hit.Detector);
    sorted_hits[detector].push_back(addressof(hit));
  }

  // apply calibration (this alters the given detectorRead!)
  for(const auto& calib : calibrations) {
    calib->ApplyTo(sorted_hits);
  }


  return nullptr;
}
