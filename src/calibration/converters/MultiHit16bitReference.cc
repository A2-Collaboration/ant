#include "MultiHit16bitReference.h"

#include "tree/TDetectorReadHit.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

const double converter::MultiHit16bitReference::CATCH_TDC_Gain = 0.1171;
const double converter::MultiHit16bitReference::V1190_TDC_Gain = 0.1;


void converter::MultiHit16bitReference::ApplyTo(const readhits_t& hits, extrahits_t&) {
    ReferenceTiming = numeric_limits<double>::quiet_NaN();

    // search for reference timing
    const auto& refhits = hits.get_item(ReferenceChannel.DetectorType);

    const auto comparer = [this] (TDetectorReadHit const * hit) {
        return hit->ChannelType == ReferenceChannel.ChannelType &&
                hit->Channel == ReferenceChannel.Channel;
    };
    const auto it_refhit = find_if(refhits.cbegin(), refhits.cend(), comparer);
    if(it_refhit == refhits.end())
        return;
    // use the same converter for the reference hit
    const vector<double> refhit_timings = ConvertWithFactorAndOffset((*it_refhit)->RawData, Gain, 0);
    // expect exactly one reference hit
    if(refhit_timings.size() != 1)
        return;
    // finally found it, remember it
    ReferenceTiming = refhit_timings.front();
}
