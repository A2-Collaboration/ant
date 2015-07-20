#include "IntegralTAPS.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDetectorRead.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

void IntegralTAPS::ProcessEvent(const Event &)
{

}

void IntegralTAPS::Finish()
{

}

void IntegralTAPS::ShowResult()
{

}

vector<double> IntegralTAPS::convert(const vector<uint8_t>& rawData) {
  if(rawData.size() != 2)
    return {};
  const uint16_t* value = reinterpret_cast<const uint16_t*>(&rawData[0]);
  return vector<double>(1, *value);
}

void IntegralTAPS::ApplyTo(const map< Detector_t::Type_t, list< TDetectorReadHit* > >& hits)
{
  // search for to be calibrated timings
  const auto it_dethits = hits.find(Detector_t::Type_t::TAPS);
  if(it_dethits == hits.end())
    return;

  const auto& dethits = it_dethits->second;

  // now calibrate the integrals (ignore any other kind of hits)
  for(TDetectorReadHit* dethit : dethits) {
    if(dethit->GetChannelType() != Channel_t::Type_t::Integral)
      continue;
    dethit->Values = convert(dethit->RawData);
    // apply offset to each of the values (might be multihit)
    const auto apply_gain = [] (double& v) {
      v *= 0.30;
      /// \todo use channel dependent gains here
    };
    for_each(dethit->Values.begin(), dethit->Values.end(),
             apply_gain);
  }

}

