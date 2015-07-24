#pragma once

#include "reconstruct/Reconstruct_traits.h"
#include "MultiHit16bit.h"

#include <limits>

namespace ant {
namespace calibration {
namespace converter {

struct CATCH_TDC : MultiHit16bit, CalibrationApply_traits {

    CATCH_TDC(const LogicalChannel_t& referenceChannel) :
        ReferenceChannel(referenceChannel),
        ReferenceTiming(std::numeric_limits<double>::quiet_NaN()),
        CATCH_to_nanoseconds(0.1171)
    {}

    virtual std::vector<double> Convert(const vector<uint8_t>& rawData) const override
    {
        // we can only convert if we have a reference hit timing
        if(std::isnan(ReferenceTiming))
            return {};
        // CATCH TDCs, the gain is fixed
        return ConvertWithFactorAndOffset(rawData, CATCH_to_nanoseconds, ReferenceTiming);
    }

    virtual void ApplyTo(const readhits_t& hits, extrahits_t&) override;

    virtual void ApplyTo(std::unique_ptr<TEvent>&) override {}

private:
    LogicalChannel_t ReferenceChannel;
    double ReferenceTiming;
    const double CATCH_to_nanoseconds;

};

}}} // namespace ant::calibration::converter