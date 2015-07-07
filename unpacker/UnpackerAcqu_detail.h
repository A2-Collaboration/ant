#ifndef UNPACKERACQU_DETAIL_H
#define UNPACKERACQU_DETAIL_H

#include <vector>
#include <memory>
#include <string>
#include <deque>
#include <ctime>
#include <cstdint>


namespace ant {

/**
 * @brief The UnpackerAcquFile class
 *
 * Base class for file access management of acqu files
 */
class RawFileReader;
class TDataRecord;
class THeaderInfo;
class UnpackerAcquConfig;

class UnpackerAcquFileFormat {
public:

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
      std::deque< std::unique_ptr<TDataRecord> >& queue
      );

  /**
   * @brief FillEvents fills the given queue with more TDataRecord items (if any left)
   * @param queue
   */
  virtual void FillEvents(std::deque< std::unique_ptr<TDataRecord> >& queue) = 0;

protected:
  virtual size_t SizeOfHeader() const = 0;
  virtual bool InspectHeader(const std::vector<uint32_t>& buffer) const = 0;
  virtual void Setup(std::unique_ptr<RawFileReader>&& reader_,
                     std::vector<std::uint32_t>&& buffer_) = 0;
  virtual void FillHeader(std::deque< std::unique_ptr<TDataRecord> >& queue) = 0;
};

// the derived file format classes
// have their own namespace
namespace unpacker {
namespace acqu {


// FileFormatBase provides a common class for Mk1/Mk2 formats (so far)
class FileFormatBase : public UnpackerAcquFileFormat {
protected:
  std::unique_ptr<RawFileReader> reader;
  std::vector<std::uint32_t> buffer;

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
    unsigned RecordLength;
  };

  Info info;
  std::uint32_t ID_upper; // upper part of UID, set by BuildTHeaderInfo
  std::unique_ptr<UnpackerAcquConfig> config;

  virtual void Setup(std::unique_ptr<RawFileReader>&& reader_,
                     std::vector<std::uint32_t>&& buffer_) override;
  virtual void FillHeader(std::deque< std::unique_ptr<TDataRecord> >& queue) override;
  virtual void FillInfo() = 0;

private:
  std::unique_ptr<THeaderInfo> BuildTHeaderInfo();
};

class FileFormatMk1 : public FileFormatBase {

  // UnpackerAcquFile interface
protected:
  virtual size_t SizeOfHeader() const override;
  virtual bool InspectHeader(const std::vector<std::uint32_t>& buffer) const override;
  virtual void FillEvents(std::deque<std::unique_ptr<TDataRecord> > &queue) override;
  virtual void FillInfo() override;
};

class FileFormatMk2 : public FileFormatBase {

  // UnpackerAcquFile interface
protected:
  virtual size_t SizeOfHeader() const override;
  virtual bool InspectHeader(const std::vector<std::uint32_t> &buffer) const override;
  virtual void FillEvents(std::deque<std::unique_ptr<TDataRecord> > &queue) override;
  virtual void FillInfo() override;

};

}} // namespace unpacker::acqu


} // namespace ant

#endif // UNPACKERACQU_DETAIL_H
