#include "TimingCATCH.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDetectorRead.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

void TimingCATCH::ProcessEvent(const Event &)
{

}

void TimingCATCH::Finish()
{

}

void TimingCATCH::ShowResult()
{

}

vector<double> convertToTiming(const vector<uint8_t>& rawData) {
  if(rawData.size() % 2 != 0)
    return {};
  vector<double> ret(rawData.size()/2);
  for(size_t i=0;i<ret.size();i++) {
    const uint16_t* rawVal = reinterpret_cast<const uint16_t*>(&rawData[2*i]);
    ret[i] = *rawVal * 0.1171; // convert it to nanoseconds
  }
  return ret;
}

void TimingCATCH::ApplyTo(const map< Detector_t::Type_t, list< TDetectorReadHit* > >& hits)
{
  // search for reference timing
  const auto it_refhits = hits.find(ReferenceChannel.DetectorType);
  if(it_refhits == hits.end())
    return;
  const list<TDetectorReadHit*>& refhits = it_refhits->second;
  auto comparer = [this] (TDetectorReadHit* hit) {
    return hit->GetChannelType() == ReferenceChannel.ChannelType &&
        hit->Channel == ReferenceChannel.Channel;
  };
  const auto it_refhit = find_if(refhits.cbegin(), refhits.cend(), comparer);
  if(it_refhit == refhits.end())
    return;

  //
  const vector<double> refhit_timings = convertToTiming((*it_refhit)->RawData);
  if(refhit_timings.size() != 1)
    return;
  const double refhit_timing = refhit_timings[0];

  // search for to be calibrated timings
  const auto it_dethits = hits.find(DetectorType);
  if(it_dethits == hits.end())
    return;

  const auto& dethits = it_dethits->second;

  // now calibrate the timings (ignore any other kind of hits)
  for(TDetectorReadHit* dethit : dethits) {
    if(dethit->GetChannelType() != Channel_t::Type_t::Timing)
      continue;
    dethit->Values = convertToTiming(dethit->RawData);
    // apply offset to each of the values (might be multihit)
    auto apply_offset = [refhit_timing] (double& v) {
      v -= refhit_timing;
      /// \todo apply channel dependent offset as well here
    };
    for_each(dethit->Values.begin(), dethit->Values.end(),
             apply_offset);
  }

}

