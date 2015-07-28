#include "ScalerFrequency.h"

#include "tree/TDetectorRead.h"
#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

void converter::ScalerFrequency::ApplyTo(const readhits_t& hits, extrahits_t&) {
    ReferenceCounts = numeric_limits<double>::quiet_NaN();

    // search for reference scaler
    const auto& refhits = hits.get_item(ReferenceScaler.DetectorType);

    const auto comparer = [this] (TDetectorReadHit const * hit) {
        return hit->GetChannelType() == ReferenceScaler.ChannelType &&
                hit->Channel == ReferenceScaler.Channel;
    };
    const auto it_refhit = find_if(refhits.cbegin(), refhits.cend(), comparer);
    if(it_refhit == refhits.end())
        return;

    ReferenceCounts = Convert32bit((*it_refhit)->RawData);
}

double converter::ScalerFrequency::Convert32bit(const vector<uint8_t>& rawData)
{
    if(rawData.size() != 4)
        return numeric_limits<double>::quiet_NaN();
    return *reinterpret_cast<const uint32_t*>(addressof(rawData[0]));
}


std::vector<double> ant::calibration::converter::ScalerFrequency::Convert(const vector<uint8_t>& rawData) const
{
    // we can only convert if we have a reference hit timing
    if(std::isnan(ReferenceCounts))
        return {};

    const double counts = Convert32bit(rawData);

    if(std::isnan(counts))
        return {};

    return { ReferenceFrequency * counts / ReferenceCounts };
}
