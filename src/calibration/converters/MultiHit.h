#pragma once

#include "Calibration.h"

#include <vector>
#include <cstdint>
#include <functional>

namespace ant {
namespace calibration {
namespace converter {

template<typename T>
struct MultiHit : Calibration::Converter {


    virtual std::vector<double> Convert(const std::vector<uint8_t>& rawData) const override
    {
        // just convert T to double
        return ConvertRaw<double>(rawData);
    }

protected:
    template<typename U = T>
    static std::vector<U> ConvertRaw(const std::vector<std::uint8_t>& rawData)
    {
        constexpr std::size_t wordsize = sizeof(T)/sizeof(std::uint8_t);
        if(rawData.size() % wordsize  != 0)
            return {};
        std::vector<U> ret(rawData.size()/wordsize);
        for(size_t i=0;i<ret.size();i++) {
            const T* rawVal = reinterpret_cast<const T*>(std::addressof(rawData[wordsize*i]));
            ret[i] = static_cast<U>(*rawVal);
        }
        return ret;
    }
};

}}} // namespace ant::calibration::converter
