#pragma once

#include "MultiHit.h"
#include "reconstruct/Reconstruct_traits.h"
#include "tree/TDetectorReadHit.h"

#include <limits>

namespace ant {
namespace calibration {
namespace converter {

struct Gains {
    const static double CATCH_TDC;
    const static double V1190_TDC;
};

template<typename T>
struct MultiHitReference : MultiHit<T>, ReconstructHook::DetectorReadHits {



    MultiHitReference(const LogicalChannel_t& referenceChannel,
                      double gain) :
        ReferenceChannel(referenceChannel),
        ReferenceHits(),
        Gain(gain)
    {}

    virtual std::vector<double> Convert(const std::vector<uint8_t>& rawData) const override
    {
        // we can only convert if we have a reference hit timing
        if(ReferenceHits.size() != 1)
            return {};
        const auto refHit = ReferenceHits.front();
        auto hits = MultiHit<T>::template ConvertRaw<double>(rawData);
        /// \todo think about hit/refHit overflow here?
        for(auto& hit : hits)
            hit = (hit - refHit)*Gain;
        return hits;
    }

    virtual void ApplyTo(const readhits_t& hits) override {
        ReferenceHits.resize(0);

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
        ReferenceHits = MultiHit<T>::template ConvertRaw<T>((*it_refhit)->RawData);
    }

protected:
    const LogicalChannel_t ReferenceChannel;
    std::vector<T> ReferenceHits; // extracted in ApplyTo
    const double Gain;
};

}}} // namespace ant::calibration::converter
