#pragma once

#include "Calibration.h"

#include <vector>
#include <cstdint>

namespace ant {
namespace calibration {
namespace converter {

struct MultiHit16bit : Calibration::Converter {


    virtual std::vector<double> Convert(const std::vector<uint8_t>& rawData) const override
    {
        return ConvertWithFactorAndOffset(rawData, 1.0, 0.0);
    }

protected:
    static std::vector<double> ConvertWithFactorAndOffset(
            const std::vector<uint8_t>& rawData,
            const double factor,
            const double offset
            ) {
        if(rawData.size() % 2 != 0)
            return {};
        std::vector<double> ret(rawData.size()/2);
        for(size_t i=0;i<ret.size();i++) {
          const uint16_t* rawVal = reinterpret_cast<const uint16_t*>(std::addressof(rawData[2*i]));
          ret[i] = *rawVal * factor - offset;
        }
        return ret;
    }
};

}}} // namespace ant::calibration::converter
