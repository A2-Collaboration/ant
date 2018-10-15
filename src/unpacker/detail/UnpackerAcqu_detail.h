#pragma once

#include "tree/TUnpackerMessage.h"
#include "UnpackerAcqu.h" // UnpackerAcquConfig

#include "base/std_ext/mapped_vectors.h"

#include <cstdint>
#include <ctime>
#include <list>
#include <memory>
#include <string>
#include <vector>

namespace ant {


/**
 * @brief The UnpackerAcquFile class
 *
 * Base class for file access management of acqu files
 */
class RawFileReader;
struct TEvent;
struct TSlowControl;

class UnpackerAcquFileFormat {
public:

    using queue_t = std::list<TEvent>; // std::list supports splice

    /**
      * @brief Get a suitable instance for the given filename
      * @param filename the file to read
      * @return the instance, or nullptr if nothing found
      *
      * Throws exception if something unusual is encountered.
      */
    static std::unique_ptr<UnpackerAcquFileFormat> Get(const std::string& filename);

    /**
      * @brief FillEvents fills the given queue with more TEvent items (if any left)
      * @param queue
      */
    virtual void FillEvents(queue_t& queue) noexcept = 0;

    virtual ~UnpackerAcquFileFormat();

    virtual double PercentDone() const =0;

protected:
    virtual size_t SizeOfHeader() const = 0;
    virtual bool InspectHeader(const std::vector<uint32_t>& buffer) const = 0;
    virtual void Setup(std::unique_ptr<RawFileReader>&& reader_,
                       std::vector<std::uint32_t>&& buffer_) = 0;
};

// the derived file format classes
// have their own namespace
namespace unpacker {
namespace acqu {

// FileFormatBase provides a common class for Mk1/Mk2 formats
class FileFormatBase : public UnpackerAcquFileFormat {
public:
    virtual ~FileFormatBase();

    virtual double PercentDone() const override;

private:
    std::unique_ptr<RawFileReader> reader;
    std::vector<std::uint32_t>     buffer;
    // messages must be buffered during event unpacking,
    // but in order to have LogMessage() const,
    // the storage must be mutable
    mutable std::vector<TUnpackerMessage>  messages;
    signed trueRecordLength;
    unsigned nUnpackedBuffers;
    unsigned nEventsInBuffer;
    time_t GetTimeStamp();
protected:

    using reader_t = decltype(reader);
    using buffer_t = decltype(buffer);
    using it_t = buffer_t::const_iterator;

    // contains what we now about the file
    struct Info {
        struct HardwareModule {
            /// \todo think about module types, such as scaler normal tdc/qdc
            std::string Identifier; // some more or less unique identifier of the module
            unsigned Index; // non-unique! Within VME it is however
            unsigned Bits;
            unsigned FirstRawChannel;
            unsigned NRawChannels;
        };
        std::vector<HardwareModule> Modules; // appear only in normal event

        std::tm Time;
        std::string Description;
        std::string RunNote;
        std::string OutFile;
        unsigned RunNumber;
        unsigned RecordLength; // Record length according to header (might not be correct, see trueRecordLength)

        enum class Format_t {
            Mk1, Mk2
        };
        Format_t Format;

    };

    Info info;

    TID id;
    unsigned AcquID_last = 0;

    using scalers_t = std::map<uint32_t, std::vector<uint32_t> >;

    // we so some more effort for the hits,
    // especially keeping storage_hits over multiple
    // events makes it considerably faster
    std::vector<UnpackerAcquConfig::hit_mapping_t> hit_mappings;
    using hit_mappings_ptr_t = std::vector< std::vector< const UnpackerAcquConfig::hit_mapping_t* > >;
    hit_mappings_ptr_t hit_mappings_ptr;
    using hit_storage_t = std_ext::mapped_vectors<std::uint16_t, std::uint16_t>;
    hit_storage_t hit_storage;

    using scaler_mappings_t = std::vector<UnpackerAcquConfig::scaler_mapping_t>;
    scaler_mappings_t scaler_mappings;


    // this class already implements some stuff
    void Setup(reader_t&& reader_, buffer_t&& buffer_) override;
    void FillEvents(queue_t& queue) noexcept override;

    // unpacker messages handling
    void LogMessage(TUnpackerMessage::Level_t level,
                    const std::string& msg, bool emit_warning = false) const;
    void AppendMessagesToEvent(TEvent& event) const;

    // Mk1/Mk2 specific methods
    virtual void FillInfo(reader_t& reader, buffer_t& buffer, Info& info) = 0;
    virtual void FillFirstDataBuffer(reader_t& reader, buffer_t& buffer) const = 0;
    virtual void UnpackEvent(TEventData& eventdata, it_t& it, const it_t& it_endbuffer, bool& good) noexcept = 0;

    // things shared by Mk1/Mk2
    bool UnpackDataBuffer(queue_t& queue, it_t& it, const it_t& it_endbuffer) noexcept;

    std::uint32_t GetDataBufferMarker() const;
    bool SearchFirstDataBuffer(reader_t& reader, buffer_t& buffer, size_t offset) const;
    bool FindFirstDataBuffer(reader_t& reader, buffer_t& buffer,
                             const size_t max_multiplier = 32,
                             const bool assert_multiplicity = true) const;
    static void FillDetectorReadHits(const hit_storage_t& hit_storage, const hit_mappings_ptr_t& hit_mappings_ptr,
                                     std::vector<TDetectorReadHit>& hits) noexcept;
    static void FillSlowControls(const scalers_t& scalers, const scaler_mappings_t& scaler_mappings,
                                 std::vector<TSlowControl>& slowcontrols) noexcept;

};

}} // namespace unpacker::acqu


} // namespace ant

