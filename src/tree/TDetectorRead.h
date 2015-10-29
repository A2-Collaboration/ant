#pragma once

#include "TDataRecord.h"

#ifndef __CINT__
#include "base/Detector_t.h"
#include <iomanip>
#include <sstream>
#endif

namespace ant {

#ifndef __CINT__
struct TDetectorReadHit  : printable_traits
#else
struct TDetectorReadHit
#endif
{
    std::uint8_t  DetectorType;
    std::uint8_t  ChannelType;
    std::uint32_t Channel;

    std::vector<std::uint8_t>   RawData;

    std::vector<double> Values;
    std::vector<bool>   ValueBits;

    const char* GetDetectorAsString() const;
    const char* GetTypeAsString() const;

#ifndef __CINT__
    TDetectorReadHit(const LogicalChannel_t& element,
                     const std::vector<std::uint8_t>& rawData) :
        DetectorType(static_cast<std::uint8_t>(element.DetectorType)),
        ChannelType(static_cast<std::uint8_t>(element.ChannelType)),
        Channel(element.Channel),
        RawData(rawData),
        Values(),
        ValueBits()
    {
        static_assert(sizeof(Channel)>=sizeof(element.Channel),
                      "LogicalElement_t::Channel does not fit into TDetecorReadHit::Channel");
    }

    TDetectorReadHit(const LogicalChannel_t& element,
                     const std::vector<double>& values) :
        DetectorType(static_cast<std::uint8_t>(element.DetectorType)),
        ChannelType(static_cast<std::uint8_t>(element.ChannelType)),
        Channel(element.Channel),
        RawData(),
        Values(values),
        ValueBits()
    {
        static_assert(sizeof(Channel)>=sizeof(element.Channel),
                      "LogicalElement_t::Channel does not fit into TDetecorReadHit::Channel");
    }

    Channel_t::Type_t GetChannelType() const {
        return static_cast<Channel_t::Type_t>(ChannelType);
    }

    Detector_t::Type_t GetDetectorType() const {
        return static_cast<Detector_t::Type_t>(DetectorType);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        s << "Hit Detector="
             //<< std::left << std::setfill(' ') << std::setw(4)
          << GetDetectorAsString() << std::right
          << " Channel="
          << std::setw(3)
          << Channel
          << " Type="
          << GetTypeAsString();

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
#ifndef ANT_TREE_ENABLE_COPY_CTOR
    TDetectorReadHit(const TDetectorReadHit&) = delete;
    TDetectorReadHit& operator=(const TDetectorReadHit&) = delete;
    TDetectorReadHit(TDetectorReadHit&&) = default;
    TDetectorReadHit& operator=(TDetectorReadHit&&) = default;
#endif
#endif

    TDetectorReadHit() :
        DetectorType(),
        ChannelType(),
        Channel(),
        RawData(),
        Values(),
        ValueBits() {}
    virtual ~TDetectorReadHit() {}
    ClassDef(TDetectorReadHit, ANT_UNPACKER_ROOT_VERSION)
};

struct TDetectorRead : TDataRecord
{
    TDetectorRead(const TID& id) :
        TDataRecord(id),
        Hits()
    {}

    std::vector<TDetectorReadHit> Hits;

    virtual void Clear() {
        Hits.resize(0);
    }

#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TDetectorRead ID=" << ID << " Hits=" << Hits.size() << '\n';
        for(size_t i=0;i<Hits.size();i++) {
            s << "  i="
              << std::setw(3)
              << i << " "
              << Hits[i] << '\n';
        }
        return s;
    }
#ifndef ANT_TREE_ENABLE_COPY_CTOR
    TDetectorRead(const TDetectorRead&) = delete;
    TDetectorRead& operator=(const TDetectorRead&) = delete;
#endif
#endif

    TDetectorRead() : TDataRecord(), Hits() {}
    ClassDef(TDetectorRead, ANT_UNPACKER_ROOT_VERSION)

};

#ifndef __CINT__
inline const char* TDetectorReadHit::GetDetectorAsString() const {
    return Detector_t::ToString(GetDetectorType());
}
inline const char* TDetectorReadHit::GetTypeAsString() const {
    return Channel_t::ToString(GetChannelType());
}
#endif

}
