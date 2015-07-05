#ifndef UNPACKERACQU_DETAIL_H
#define UNPACKERACQU_DETAIL_H

#include "printable.h"
#include "stl_helpers.h"
#include <memory>
#include <string>
#include <deque>

namespace ant {

/**
 * @brief The UnpackerAcquFile class
 *
 * Base class for file access management of acqu files
 */

class RawFileReader;
class TDataRecord;



class UnpackerAcquFileFormat {
public:

  /**
   * @brief Get
   * @param filename
   * @param queue
   * @return
   *
   * This factory method returns a fully setup reader. It might already fill
   * the quere with some events about the header or some unpacker messages, for example.
   */
  static std::unique_ptr<UnpackerAcquFileFormat> Get(const std::string& filename,
                                                     std::deque< std::unique_ptr<TDataRecord> >& queue
                                                     );

  virtual void FillEvents(std::deque< std::unique_ptr<TDataRecord> >& queue) = 0;

protected:

  virtual size_t SizeOfHeader() const = 0;
  virtual bool InspectHeader(const std::vector<uint32_t>& buffer) const = 0;
  virtual void Setup(std::unique_ptr<RawFileReader>&& reader_,
                     std::vector<uint32_t>&& buffer_) = 0;
  virtual void FillHeader(std::deque< std::unique_ptr<TDataRecord> >& queue) = 0;
};

// the derived file format classes
// have their own namespace
namespace unpacker {
namespace acqu {


class FileFormatBase : public UnpackerAcquFileFormat {
protected:
  std::unique_ptr<RawFileReader> reader;
  std::vector<uint32_t> buffer;

  virtual void Setup(std::unique_ptr<RawFileReader>&& reader_,
                     std::vector<uint32_t>&& buffer_) override;
};

class FileFormatMk1 : public FileFormatBase {

  // UnpackerAcquFile interface
protected:
  virtual size_t SizeOfHeader() const override;
  virtual bool InspectHeader(const std::vector<uint32_t> &buffer) const override;
  virtual void FillEvents(std::deque<std::unique_ptr<TDataRecord> > &queue) override;
  virtual void FillHeader(std::deque< std::unique_ptr<TDataRecord> >& queue) override;
};

class FileFormatMk2 : public FileFormatBase {

  // UnpackerAcquFile interface
protected:
  virtual size_t SizeOfHeader() const override;
  virtual bool InspectHeader(const std::vector<uint32_t> &buffer) const override;
  virtual void FillEvents(std::deque<std::unique_ptr<TDataRecord> > &queue) override;
  virtual void FillHeader(std::deque< std::unique_ptr<TDataRecord> >& queue) override;
};

}} // namespace unpacker::acqu


} // namespace ant

#endif // UNPACKERACQU_DETAIL_H
