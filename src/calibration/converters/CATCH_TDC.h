#pragma once

#include "reconstruct/Reconstruct_traits.h"
#include "MultiHitReference.h"

#include <limits>

namespace ant {
namespace calibration {
namespace converter {

struct CATCH_TDC : MultiHitReference<std::uint16_t> {

    CATCH_TDC(const LogicalChannel_t& referenceChannel) :
        MultiHitReference(referenceChannel, Gains::CATCH_TDC)
    {}

    virtual std::vector<double> Convert(const std::vector<uint8_t>& rawData) const override
    {
        // we can only convert if we have exactly one reference hit timing
        if(ReferenceHits.size() != 1)
            return {};
        const std::int32_t refHit = ReferenceHits.front();
        // reject conversion if refhit is invalid (0xffff)
        constexpr std::uint16_t max_u16bit = std::numeric_limits<std::uint16_t>::max();
        if(refHit == max_u16bit)
            return {};

        auto rawHits = MultiHit<std::uint16_t>::template ConvertRaw<std::int32_t>(rawData);

        // the magic value was originally 62054, but
        // investigating the output of the CATCH TDC showed that 62121 seems more
        // like the "true" overflow value of the F1 chip
        constexpr std::int32_t CATCH_Overflow = 62054;

        std::vector<double> hits;
        for(const std::uint16_t rawHit : rawHits) {
            // reject invalid rawhits
            if(rawHit == max_u16bit) {
                continue;
            }
            // use the same check for overflows as Acqu does
            // it is important that rawHit/refHit are int32
            // to prevent wrap-arounds in subtraction
            auto value = static_cast<std::int32_t>(rawHit) - static_cast<std::int32_t>(refHit);
            const auto value_p = value + CATCH_Overflow;
            const auto value_m = value - CATCH_Overflow;
            value = abs(value) < abs(value_p) ? value : value_p;
            value = abs(value) < abs(value_m) ? value : value_m;
            hits.push_back(value*Gain);
        }

        return hits;
    }
};

}}} // namespace ant::calibration::converter
