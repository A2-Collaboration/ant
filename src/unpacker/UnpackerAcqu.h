#pragma once

#include "Unpacker.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/Detector_t.h"

#include <memory>
#include <list>
#include <map>
#include <vector>
#include <cstdint>
#include <limits>

namespace ant {

class UnpackerAcquFileFormat; // implemented in detail/UnpackerAcqu.h

class UnpackerAcqu : public Unpacker::Module
{
public:
    UnpackerAcqu();
    virtual bool OpenFile(const std::string& filename) override;
    virtual std::shared_ptr<TDataRecord> NextItem() noexcept override;

    class Exception : public Unpacker::Exception {
        using Unpacker::Exception::Exception; // use base class constructor
    };

private:
    std::list< std::unique_ptr<TDataRecord> > queue; // std::list supports splice
    std::unique_ptr<UnpackerAcquFileFormat> file;

};

// we define some methods here which
// the configs are required to implement
class UnpackerAcquConfig : public ExpConfig::Unpacker<UnpackerAcquConfig> {
public:

    template<typename T>
    struct RawChannel_t {
        static_assert(std::is_unsigned<T>::value, "T must be unsigned");
        T RawChannel;
        T Mask;
        static constexpr T NoMask = std::numeric_limits<T>::max();
        // provide some handy constructors (implemented below)
        RawChannel_t(const std::initializer_list<T>& l);
        RawChannel_t(const T& ch);
    };

    template<typename T>
    struct mapping_t {
        LogicalChannel_t LogicalChannel;
        std::vector< RawChannel_t<T> > RawChannels;
        // provide constructors for single raw channel mapping
        mapping_t(
                const LogicalChannel_t& logicalChannel,
                T rawChannel
                ) :
            LogicalChannel(logicalChannel),
            RawChannels{rawChannel}
        {}
        mapping_t(
                Detector_t::Type_t detector,
                Channel_t::Type_t channeltype,
                unsigned channel,
                T rawChannel
                ) :
            LogicalChannel{detector, channeltype, channel},
            RawChannels{rawChannel}
        {}
    };

    // hits (or Acqu ADC values) are always part of TDetectorRead
    struct hit_mapping_t : mapping_t<std::uint16_t> {
        using mapping_t::mapping_t; // use constructors from base class
    };

    // scalers in acqu can be handled as additional information
    // for a logical detector channel, or as TSlowControl items
    struct scaler_mapping_t : mapping_t<std::uint32_t> {
        // if non-empty, scaler is converted to TSlowControl item
        std::string SlowControlName;

        // provide constructors
        scaler_mapping_t(
                const std::string& slowControlName,
                Detector_t::Type_t detector,
                unsigned channel,
                std::uint32_t rawChannel
                ) :
            mapping_t(detector, Channel_t::Type_t::Scaler, channel, rawChannel),
            SlowControlName(slowControlName)
        {
            if(SlowControlName.empty())
                throw std::runtime_error("Do not use this ctor with empty slowControlName");
        }
        scaler_mapping_t(
                Detector_t::Type_t detector,
                unsigned channel,
                std::uint32_t rawChannel
                ) :
            mapping_t(detector, Channel_t::Type_t::Scaler, channel, rawChannel),
            SlowControlName()
        {}
    };

    virtual void BuildMappings(
            std::vector<hit_mapping_t>& hit_mappings,
            std::vector<scaler_mapping_t>& scaler_mappings
            ) const = 0;
};

// define the templated constructors here to keep the class definition clean
template<typename T>
inline UnpackerAcquConfig::RawChannel_t<T>::RawChannel_t(const std::initializer_list<T> &l) {
    if(l.size()==2) {
        const std::vector<T> v(l);
        RawChannel = v[0];
        Mask = v[1];
    }
    else
        throw std::runtime_error("RawChannel_t can only be initialized with 2 values.");
}

template<typename T>
inline UnpackerAcquConfig::RawChannel_t<T>::RawChannel_t(const T &ch)
{
    RawChannel = ch;
    Mask = RawChannel_t<T>::NoMask;
}


} // namespace ant
