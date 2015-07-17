#include "IntegralSADC.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDetectorRead.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

void IntegralSADC::ProcessEvent(const Event &)
{

}

void IntegralSADC::Finish()
{

}

void IntegralSADC::ShowResult()
{

}

vector<double> convertToValue(const vector<uint8_t>& rawData) {
  if(rawData.size() != 6) // need three 16bit values
    return {};

  vector<double> ret(1);
  const uint16_t* pedestal = reinterpret_cast<const uint16_t*>(&rawData[0]);
  const uint16_t* signal = reinterpret_cast<const uint16_t*>(&rawData[2]);

  return vector<double>(1, *signal - *pedestal);
}

void IntegralSADC::ApplyTo(const map< Detector_t::Type_t, list< TDetectorReadHit* > >& hits)
{
  // search for to be calibrated timings
  const auto it_dethits = hits.find(DetectorType);
  if(it_dethits == hits.end())
    return;

  const auto& dethits = it_dethits->second;

  // now calibrate the integrals (ignore any other kind of hits)
  for(TDetectorReadHit* dethit : dethits) {
    if(dethit->GetChannelType() != Channel_t::Type_t::Integral)
      continue;
    dethit->Values = convertToValue(dethit->RawData);
    // apply offset to each of the values (might be multihit)
    const auto apply_gain = [] (double& v) {
      v *= 0.07;
      /// \todo use channel dependent gains here
    };
    for_each(dethit->Values.begin(), dethit->Values.end(),
             apply_gain);
  }

}

