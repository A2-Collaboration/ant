#include "CATCH_TDC.h"

#include "tree/TDetectorRead.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

void converter::CATCH_TDC::ApplyTo(TDetectorRead&, const readhits_t& hits) {
    ReferenceTiming = numeric_limits<double>::quiet_NaN();
    // search for reference timing
    const auto it_refhits = hits.find(ReferenceChannel.DetectorType);
    if(it_refhits == hits.end())
        return;
    const auto& refhits = it_refhits->second;
    const auto comparer = [this] (TDetectorReadHit const * hit) {
        return hit->GetChannelType() == ReferenceChannel.ChannelType &&
                hit->Channel == ReferenceChannel.Channel;
    };
    const auto it_refhit = find_if(refhits.cbegin(), refhits.cend(), comparer);
    if(it_refhit == refhits.end())
        return;
    // use the same converter for the reference hit
    const vector<double> refhit_timings = ConvertWithFactorAndOffset((*it_refhit)->RawData, CATCH_to_nanoseconds, 0);
    // expect exactly one reference hit
    if(refhit_timings.size() != 1)
        return;
    // finally found it, remember it
    ReferenceTiming = refhit_timings[0];
}
