#pragma once

#include "calibration/Calibration.h"

#include <vector>
#include <cstdint>

namespace ant {
namespace calibration {
namespace converter {

struct GeSiCa_SADC : Calibration::Converter {


    virtual std::vector<double> Convert(const vector<uint8_t>& rawData) const override
    {
        if(rawData.size() != 6) // expect three 16bit values
          return {};

        const double pedestal = *reinterpret_cast<const uint16_t*>(&rawData[0]);
        const double signal = *reinterpret_cast<const uint16_t*>(&rawData[2]);

        // return vector with size 1 and pedestal subtracted signal
        return {signal - pedestal};
    }
};

}}} // namespace ant::calibration::converter
