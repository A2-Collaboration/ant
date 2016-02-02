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

    std::vector<std::uint8_t>   RawData;

    std::vector<double> Converted;
    std::vector<double> Values;
    std::vector<bool>   ValueBits;

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

    TDetectorReadHit(const LogicalChannel_t& element,
                     const std::vector<double>& values) :
        DetectorType(element.DetectorType),
        ChannelType(element.ChannelType),
        Channel(element.Channel),
        RawData(),
        Values(values),
        ValueBits()
    {
    }

    TDetectorReadHit() {}
    virtual ~TDetectorReadHit() {}

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
            std::ostringstream s_values;
            for(double c : Values) {
                s_values << c << " ";
            }
            s << " Values=" << s_values.str();
        }

        return s;
    }

    TDetectorReadHit(const TDetectorReadHit&) = delete;
    TDetectorReadHit& operator=(const TDetectorReadHit&) = delete;
    TDetectorReadHit(TDetectorReadHit&&) = default;
    TDetectorReadHit& operator=(TDetectorReadHit&&) = default;


};

}
