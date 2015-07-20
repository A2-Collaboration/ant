#include "TAPS.h"

#include "detail/TAPS_2013_BaF2_elements.h"
#include "detail/TAPS_2013_PbWO4_elements.h"


#include "tree/THeaderInfo.h"
#include "base/std_ext.h"

#include <iostream>
#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;


void TAPS::BuildMappings(
    vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
    vector<UnpackerAcquConfig::scaler_mapping_t> &) const
{
  const auto& baf2s = GetBaF2Elements();
  const auto& pbwo4s = GetPbWO4Elements();


}

void TAPS::BuildClusterElements()
{
  // get the elements from the derived class
  // we cannout use GetClusterElement somehow because
  // it returns a const ptr
  auto& baf2s = GetBaF2Elements();
  auto& pbwo4s = GetPbWO4Elements();

  clusterelements.resize(baf2s.size()+pbwo4s.size(), nullptr);

  // apply the z-position depending on Cherenkov
  // and build the clusterelements
  // we assume that the channel elements are consecutive
  const double zpos = CherenkovInstalled ? 174.2 : 145.7;
  for(auto& baf2 : baf2s) {
    baf2.Position.SetZ(zpos);
    clusterelements[baf2.Channel] = addressof(baf2);
  }
  for(auto& pbwo4 : pbwo4s) {
    pbwo4.Position.SetZ(zpos);
    clusterelements[pbwo4.Channel] = addressof(pbwo4);
  }

  // check that every element was actually set
  for(auto elem : clusterelements) {
    assert(elem != nullptr);
  }
}

bool TAPS_2013::Matches(const THeaderInfo& headerInfo) const {
  return std_ext::time_after(headerInfo.Timestamp, "2013-11-01");
}


