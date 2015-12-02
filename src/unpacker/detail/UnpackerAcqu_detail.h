#pragma once

#include "tree/TUnpackerMessage.h" // TUnpackerMessage::Level_t
#include "UnpackerAcqu.h" // UnpackerAcquConfig

#include "base/mapped_vectors.h"

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
struct TDataRecord;
struct THeaderInfo;

class UnpackerAcquFileFormat {
public:

    using queue_t = std::list< std::unique_ptr<TDataRecord> >;

    /**
   * @brief Get a suitable instance for the given filename
   * @param filename the file to read
   * @param queue for possible ant::T* messages during setup
   * @return the instance, or nullptr if nothing found
   *
   * This factory method returns a fully setup reader. It might already fill
   * the queue with some events about the header or some unpacker messages, for example.
   *
   * Throws exception if something unusual is encountered.
   *
   */
    static std::unique_ptr<UnpackerAcquFileFormat> Get(
            const std::string& filename,
            queue_t& queue
            );

    /**
   * @brief FillEvents fills the given queue with more TDataRecord items (if any left)
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
    virtual void FillHeader(queue_t& queue) = 0;
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
    std::vector<std::uint32_t> buffer;
    signed trueRecordLength;
    unsigned unpackedBuffers;
    std::unique_ptr<THeaderInfo> BuildTHeaderInfo();
    time_t GetTimeStamp();
protected:

    using reader_t = decltype(reader);
    using buffer_t = decltype(buffer);
    using it_t = buffer_t::const_iterator;

    // helper class to finally create the THeaderInfo items
    // the fields are still very similar to the Acqu header info fields
    // but the knowledge what data format it actually was created from
    // is already gone, yay
    struct Info {
        struct HardwareModule {
            std::string Identifier; // some more or less unique identifier of the module
            unsigned Index; // non-unique! Within VME it is however
            unsigned Bits;
            unsigned FirstRawChannel;
            unsigned NRawChannels;
        };
        std::vector<HardwareModule> ADCModules; // appear only in normal event
        std::vector<HardwareModule> ScalerModules; // appear only in scaler blocks

        std::tm Time;
        std::string Description;
        std::string RunNote;
        std::string OutFile;
        unsigned RunNumber;
        unsigned RecordLength; // Record length according to header (might not be correct, see trueRecordLength)
    };

    Info info;
    TID id;
    unsigned AcquID_last = 0;

    // we so some more effort for the hits,
    // especially keeping storage_hits over multiple
    // events makes it considerably faster
    std::vector<UnpackerAcquConfig::hit_mapping_t> hit_mappings;
    std::vector< std::vector< const UnpackerAcquConfig::hit_mapping_t* > > hit_mappings_ptr;
    using hits_t = std_ext::mapped_vectors<uint16_t, uint16_t>;
    hits_t hit_storage;

    std::vector<UnpackerAcquConfig::scaler_mapping_t> scaler_mappings;


    // this class already implements some stuff
    void Setup(reader_t&& reader_, buffer_t&& buffer_) override;
    void FillHeader(queue_t& queue) override;
    void FillEvents(queue_t& queue) noexcept override;
    void LogMessage(queue_t& queue,
                    TUnpackerMessage::Level_t level,
                    const std::string &msg) const;

    // Mk1/Mk2 specific methods
    virtual void FillInfo(reader_t& reader, buffer_t& buffer, Info& info) const = 0;
    virtual void FillFirstDataBuffer(queue_t& queue, reader_t& reader, buffer_t& buffer) const = 0;
    virtual bool UnpackDataBuffer(queue_t &queue, it_t& it, const it_t& it_endbuffer) noexcept = 0;

};

}} // namespace unpacker::acqu


} // namespace ant

