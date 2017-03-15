#pragma once

#include "base/Detector_t.h"
#include <iomanip>
#include <sstream>

namespace ant {

struct TDetectorReadHit  : printable_traits
{
    Detector_t::Type_t DetectorType;
    Channel_t::Type_t  ChannelType;
    std::uint32_t      Channel;

    // represents some arbitrary binary blob
    std::vector<std::uint8_t> RawData;

    // encapsulates the possible outcomes of conversion
    // from RawData, including intermediate results (typically before calibration)
    struct Value_t : printable_traits {
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

        virtual std::ostream& Print( std::ostream& s) const override {
            if(!std::isnan(Calibrated))
                s << Calibrated;
            else if(!std::isnan(Uncalibrated))
                s << "U=" << Uncalibrated;
            else
                s << "undef";
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
    virtual ~TDetectorReadHit() = default;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(DetectorType, ChannelType, Channel, RawData, Values, ValueBits);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        s << "Hit Detector="
          << Detector_t::ToString(DetectorType) << std::right
          << " Channel="
          << std::setw(3)
          << Channel
          << " Type="
          << Channel_t::ToString(ChannelType);

        if(!RawData.empty()) {
            std::ostringstream s_rawdata;
            s_rawdata << std::hex << std::uppercase << std::setfill('0');
            for(auto i = RawData.rbegin(); i != RawData.rend(); i++) {
                s_rawdata << std::setw(2) << static_cast<int>(*i);
            }
            s << " RawData=0x" << s_rawdata.str();
        }
        if(!Values.empty()) {
            s << " Values=";
            for(auto& c : Values) {
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
