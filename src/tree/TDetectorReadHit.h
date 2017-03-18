#pragma once

#include "base/Detector_t.h"
#include <iomanip>
#include <sstream>

namespace ant {

struct TDetectorReadHit
{
    Detector_t::Type_t DetectorType;
    Channel_t::Type_t  ChannelType;
    std::uint32_t      Channel;

    // represents some arbitrary binary blob
    std::vector<std::uint8_t> RawData;

    // encapsulates the possible outcomes of conversion
    // from RawData, including intermediate results (typically before calibration)
    struct Value_t {
        explicit Value_t(double uncalibrated = std_ext::NaN)
            : Uncalibrated(uncalibrated),
              Calibrated(uncalibrated)
        {}
        double Uncalibrated; // converted value, for example useful for pedestals
        double   Calibrated; // final value, commonly used as timings or energy

        template<class Archive>
        void serialize(Archive& archive) {
            archive(Uncalibrated, Calibrated);
        }

        friend std::ostream& operator<<( std::ostream& s, const Value_t& o) {
            if(!std::isnan(o.Calibrated))
                s << o.Calibrated;
            else if(!std::isnan(o.Uncalibrated))
                s << "U=" << o.Uncalibrated;
            else
                s << "UNDEFINED";
            return s;
        }
    };

    std::vector<Value_t> Values;
    std::vector<bool>    ValueBits;

    // RawData ctor
    TDetectorReadHit(const LogicalChannel_t& element,
                     const std::vector<std::uint8_t>& rawData) :
        DetectorType(element.DetectorType),
        ChannelType(element.ChannelType),
        Channel(element.Channel),
        RawData(rawData),
        Values(),
        ValueBits()
    {
    }

    // Single (typically uncalibrated) value ctor
    TDetectorReadHit(const LogicalChannel_t& element,
                     const Value_t& value) :
        DetectorType(element.DetectorType),
        ChannelType(element.ChannelType),
        Channel(element.Channel),
        RawData(),
        Values{value},
        ValueBits()
    {
    }

    TDetectorReadHit() = default;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(DetectorType, ChannelType, Channel, RawData, Values, ValueBits);
    }

    friend std::ostream& operator<<( std::ostream& s, const TDetectorReadHit& o) {
        s << "Hit Detector="
          << Detector_t::ToString(o.DetectorType) << std::right
          << " Channel="
          << std::setw(3)
          << o.Channel
          << " Type="
          << Channel_t::ToString(o.ChannelType);

        if(!o.RawData.empty()) {
            std::ostringstream s_rawdata;
            s_rawdata << std::hex << std::uppercase << std::setfill('0');
            for(auto i = o.RawData.rbegin(); i != o.RawData.rend(); i++) {
                s_rawdata << std::setw(2) << static_cast<int>(*i);
            }
            s << " RawData=0x" << s_rawdata.str();
        }
        if(!o.Values.empty()) {
            s << " Values=";
            for(auto& c : o.Values) {
                s << c << " ";
            }
        }

        return s;
    }

    TDetectorReadHit(const TDetectorReadHit&) = delete;
    TDetectorReadHit& operator=(const TDetectorReadHit&) = delete;
    TDetectorReadHit(TDetectorReadHit&&) = default;
    TDetectorReadHit& operator=(TDetectorReadHit&&) = default;


};

}
