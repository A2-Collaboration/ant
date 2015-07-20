#include "TimingTAPS.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDetectorRead.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

void TimingTAPS::ProcessEvent(const Event &)
{

}

void TimingTAPS::Finish()
{

}

void TimingTAPS::ShowResult()
{

}

vector<double> TimingTAPS::convert(const vector<uint8_t>& rawData) {
  if(rawData.size() % 2 != 0)
    return {};
  vector<double> ret(rawData.size()/2);
  for(size_t i=0;i<ret.size();i++) {
    const uint16_t* rawVal = reinterpret_cast<const uint16_t*>(&rawData[2*i]);
    ret[i] = *rawVal;
  }
  return ret;
}

void TimingTAPS::ApplyTo(const map< Detector_t::Type_t, list< TDetectorReadHit* > >& hits)
{
  // search for to be calibrated timings
  const auto it_dethits = hits.find(Detector_t::Type_t::TAPS);
  if(it_dethits == hits.end())
    return;

  const auto& dethits = it_dethits->second;

  // now calibrate the timings (ignore any other kind of hits)
  for(TDetectorReadHit* dethit : dethits) {
    if(dethit->GetChannelType() != Channel_t::Type_t::Timing)
      continue;
    dethit->Values = convert(dethit->RawData);

    const auto apply_offset = [] (double& v) {
      v *= -0.100;
      /// \todo apply channel dependent timing correction here
    };
    for_each(dethit->Values.begin(), dethit->Values.end(),
             apply_offset);
  }

}

