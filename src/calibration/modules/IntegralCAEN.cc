#include "IntegralCAEN.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDetectorRead.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

void IntegralCAEN::ProcessEvent(const Event &)
{

}

void IntegralCAEN::Finish()
{

}

void IntegralCAEN::ShowResult()
{

}

vector<double> IntegralCAEN::convert(const vector<uint8_t>& rawData) {
    if(rawData.size() % 2 != 0)
        return {};
    vector<double> ret(rawData.size()/2);
    for(size_t i=0;i<ret.size();i++) {
        const uint16_t* rawVal = reinterpret_cast<const uint16_t*>(&rawData[2*i]);
        ret[i] = *rawVal;
    }
    return ret;
}

void IntegralCAEN::ApplyTo(const map< Detector_t::Type_t, list< TDetectorReadHit* > >& hits)
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
        dethit->Values = convert(dethit->RawData);
        // apply offset to each of the values (might be multihit)
        const auto apply_gain = [] (double& v) {
            v *= 0.014;
            /// \todo use channel dependent gains here
        };
        for_each(dethit->Values.begin(), dethit->Values.end(),
                 apply_gain);
    }

}

