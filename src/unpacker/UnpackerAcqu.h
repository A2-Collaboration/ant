#pragma once

#include "Unpacker.h"

#include "expconfig/ExpConfig.h"
#include "base/Detector_t.h"

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
    virtual TEvent NextEvent() override;
    virtual bool ProvidesSlowControl() const override { return true; }

    class Exception : public Unpacker::Exception {
        using Unpacker::Exception::Exception; // use base class constructor
    };

    virtual ~UnpackerAcqu();

    virtual double PercentDone() const override;

private:
    std::list<TEvent> queue; // std::list supports splice
    std::unique_ptr<UnpackerAcquFileFormat> file;

};

// we define some methods here which
// the configs are required to implement
class UnpackerAcquConfig {
public:

    template<typename T>
    struct RawChannel_t {
        static_assert(std::is_unsigned<T>::value, "T must be unsigned");
        T RawChannel;
        T Mask;
        static constexpr T NoMask() { return std::numeric_limits<T>::max(); }
        // provide some handy constructors
        RawChannel_t(const T& ch, const T& mask) : RawChannel(ch), Mask(mask) {}
        RawChannel_t(const T& ch) : RawChannel_t(ch, NoMask()) {}
    };

    // hits (or Acqu ADC values) are always part of TDetectorRead
    //
    struct hit_mapping_t {
        LogicalChannel_t LogicalChannel;
        std::vector< RawChannel_t< std::uint16_t > > RawChannels;
        // provide constructors for single raw channel mapping
        // (most commonly used)
        hit_mapping_t(
                const LogicalChannel_t& logicalChannel,
                std::uint16_t rawChannel
                ) :
            LogicalChannel(logicalChannel),
            RawChannels{rawChannel}
        {}
        hit_mapping_t(
                Detector_t::Type_t detector,
                Channel_t::Type_t channeltype,
                unsigned channel,
                std::uint16_t rawChannel
                ) :
            LogicalChannel{detector, channeltype, channel},
            RawChannels{rawChannel}
        {}
    };

    // scalers in acqu are always TSlowControl items
    struct scaler_mapping_t  {
        std::string SlowControlName; // must not be empty!

        // slowcontrol may have more than one
        // data entry (for example tagger scaler for each channel)
        struct entry_t {
            unsigned LogicalChannel;
            RawChannel_t< std::uint32_t > RawChannel;

            entry_t(unsigned channel, std::uint32_t rawChannel) :
                LogicalChannel(channel), RawChannel(rawChannel)
            {}
        };
        std::vector< entry_t > Entries;

        // provide constructors for simple creation
        scaler_mapping_t(
                const std::string& slowControlName,
                const std::vector< entry_t >& entries
                ) :
            SlowControlName(slowControlName),
            Entries(entries)
        {
            if(SlowControlName.empty())
                throw std::runtime_error("Empty SlowControlName provided");
        }

        scaler_mapping_t(
                const std::string& slowControlName,
                std::uint32_t rawChannel
                ) :
            scaler_mapping_t(slowControlName, {entry_t(0, rawChannel)})
        {}
    };

    virtual void BuildMappings(
            std::vector<hit_mapping_t>& hit_mappings,
            std::vector<scaler_mapping_t>& scaler_mappings
            ) const = 0;

protected:
    ~UnpackerAcquConfig() = default;
};

} // namespace ant
