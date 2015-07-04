#ifndef UNPACKERACQU_DETAIL_H
#define UNPACKERACQU_DETAIL_H

#include "printable.h"
#include "stl_helpers.h"
#include <memory>
#include <string>


namespace ant {

/**
 * @brief The UnpackerAcquFile class
 *
 * Base class for file access management of acqu files
 */

class RawFileReader;

class UnpackerAcquFileFormat {
public:
  virtual ~UnpackerAcquFileFormat() = default;
  // factory method to get concrete implementation
  static std::unique_ptr<UnpackerAcquFileFormat> Get(const std::string& filename);

  class Info {

  };
  virtual Info GetInfo() { return Info(); }

protected:
  virtual size_t SizeOfHeader() const = 0;
  virtual bool InspectHeader(const std::vector<uint32_t>& buffer) const = 0;
  virtual void Setup(std::unique_ptr<RawFileReader>&& reader, std::vector<uint32_t>&& buffer) = 0;
};

// the derived file format classes
// have their own namespace
namespace unpacker {
namespace acqu {

class FileFormatMk1 : public UnpackerAcquFileFormat {

  // UnpackerAcquFile interface
protected:
  virtual size_t SizeOfHeader() const override;
  virtual bool InspectHeader(const std::vector<uint32_t> &buffer) const override;
  virtual void Setup(std::unique_ptr<RawFileReader>&& reader, std::vector<uint32_t>&& buffer) override;
};

class FileFormatMk2 : public UnpackerAcquFileFormat {

  // UnpackerAcquFile interface
protected:
  virtual size_t SizeOfHeader() const override;
  virtual bool InspectHeader(const std::vector<uint32_t> &buffer) const override;
  virtual void Setup(std::unique_ptr<RawFileReader>&& reader, std::vector<uint32_t>&& buffer) override;
};

}} // namespace unpacker::acqu


} // namespace ant

#endif // UNPACKERACQU_DETAIL_H
