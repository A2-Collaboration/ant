#pragma once

#include "reconstruct/Reconstruct_traits.h"
#include "MultiHit16bit.h"

#include <limits>

namespace ant {
namespace calibration {
namespace converter {

struct MultiHit16bitReference : MultiHit16bit, ReconstructHook::DetectorReadHits {

    const static double CATCH_TDC_Gain;
    const static double V1190_TDC_Gain;

    MultiHit16bitReference(const LogicalChannel_t& referenceChannel, double gain) :
        ReferenceChannel(referenceChannel),
        ReferenceTiming(std::numeric_limits<double>::quiet_NaN()),
        Gain(gain)
    {}

    virtual std::vector<double> Convert(const std::vector<uint8_t>& rawData) const override
    {
        // we can only convert if we have a reference hit timing
        if(std::isnan(ReferenceTiming))
            return {};

        return ConvertWithFactorAndOffset(rawData, Gain, ReferenceTiming);
    }

    virtual void ApplyTo(const readhits_t& hits, extrahits_t&) override;

private:
    LogicalChannel_t ReferenceChannel;
    double ReferenceTiming;
    const double Gain;

};

}}} // namespace ant::calibration::converter
