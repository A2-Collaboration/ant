#pragma once

#include "Calibration.h"
#include "reconstruct/Reconstruct_traits.h"

#include <limits>

namespace ant {
namespace calibration {
namespace converter {

struct ScalerFrequency : Calibration::Converter, CalibrationApply_traits {

    ScalerFrequency(const LogicalChannel_t& referenceScaler) :
        ReferenceScaler(referenceScaler),
        ReferenceCounts(std::numeric_limits<double>::quiet_NaN()),
        ReferenceFrequency(1e6) // assume that all reference counters are 1MHz
    {}

    virtual std::vector<double> Convert(const vector<uint8_t>& rawData) const override;

    virtual void ApplyTo(const readhits_t& hits, extrahits_t&) override;

    virtual void ApplyTo(std::unique_ptr<TEvent>&) override {}

private:
    LogicalChannel_t ReferenceScaler;
    double ReferenceCounts;
    const double ReferenceFrequency;

    static double Convert32bit(const vector<uint8_t>& rawData);
};

}}} // namespace ant::calibration::converter
