#include "TAPS.h"

#include "tree/THeaderInfo.h"
#include "base/std_ext.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;

bool TAPS_2013::Matches(const THeaderInfo& headerInfo) const {
  return std_ext::time_after(headerInfo.Timestamp, "2013-11-01");
}

void TAPS_2013::BuildMappings(
    vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
    vector<UnpackerAcquConfig::scaler_mapping_t> &) const
{
  for(unsigned sector=0;sector<6;sector++) {

  }
}
