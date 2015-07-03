#ifndef RAWFILEREADER_H
#define RAWFILEREADER_H

#include <fstream>
#include <string>
#include <memory>
#include <cstdint>
#include <vector>

#include "stl_helpers.h"

// fix qtcreator highlighting...
typedef std::uint8_t uint8_t;

namespace ant {

/**
 * @brief The RawFileReader class
 * Super-simple wrapper for reading binary files,
 * even if they are compressed
 * Possible IO errors are propagated as exceptions
 */
class RawFileReader {

public:

  /**
   * @brief RawFileReader
   * @param filename
   * @param inbufsize
   *
   * Parameter inbufsize is ignored if non-compressed data is read
   */
  explicit RawFileReader(const std::string& filename, const size_t inbufsize = BUFSIZ)
  {
    if(std_ext::string_ends_with(filename, ".xz")) {
      p = std::unique_ptr<PlainBase>(new XZ(filename, inbufsize));
    }
    else {
      p = std::unique_ptr<PlainBase>(new PlainBase(filename));
    }
  }

  // methods/operators for testing the state
  //bool is_open() const { return file.is_open(); }
  //bool operator!() const { return file.operator!(); }
  explicit operator bool() const {
    return p->operator bool();
  }

  /**
   * @brief read n bytes to buffer s
   * @param s buffer (must be allocated by caller)
   * @param n
   */
  void read(char* s, std::streamsize n) {
    p->read(s,n);
  }

  /**
   * @brief gcount
   * @return number of bytes read
   */
  std::streamsize gcount() const {
    return p->gcount();
  }

  /**
   * @brief eof
   * @return true if last read went beyond end of file.
   *
   * Note that reading the exact number of bytes
   */
  bool eof() const {
    return p->eof();
  }

  class Exception : public std::runtime_error {
    using std::runtime_error::runtime_error; // use base class constructor
  };

private:

  /**
   * @brief The PlainBase class
   *
   * Plain and simple filereader, encapsulating an ifstream
   * and providing a base class for the more complicated decompression classes
   *
   * Only the really needed methods are exported
   */
  class PlainBase {
  public:
    explicit PlainBase(const std::string& filename)
      : file(filename.c_str(), std::ios::binary)
    {}

    virtual ~PlainBase() {} // derived classes may want to

    // methods/operators for testing the state
    virtual explicit operator bool() const {
      return !file.operator!(); // some older ifstream version don't implement "operator bool"
    }

    // read methods
    virtual void read(char* s, std::streamsize n) {
      file.read(s, n);
    }

    virtual bool eof() const {
      return file.eof();
    }

    virtual std::streamsize gcount() const {
      return file.gcount();
    }

  private:
    std::ifstream file;
  }; // class RawFileReader::Plain

  /**
   * @brief The XZ class
   * Decompresses the bytes before giving them to the user
   *
   * Adapted from
   * http://git.tukaani.org/?p=xz.git;a=blob_plain;f=doc/examples/02_decompress.c;hb=HEAD
   *
   */
  class XZ : public PlainBase {
  public:

    XZ(const std::string& filename, const size_t inbufsize);

    virtual ~XZ();

    virtual explicit operator bool() const {
      return PlainBase::operator bool() && !decompressFailed;
    }

    virtual void read(char *s, std::streamsize n);

    virtual std::streamsize gcount() const {
      return gcount_;
    }

    virtual bool eof() const {
      return eof_;
    }

  private:
    std::vector<uint8_t> inbuf;
    bool decompressFailed;
    std::streamsize gcount_;
    bool eof_;


#ifndef RAWFILEREADER_H_IMPL
    struct lzma_stream;
#endif
    lzma_stream* strm;
    void init_decoder();
    void cleanup();

  }; // class RawFileReader::XZ


  // private stuff for RawFileReader
  std::unique_ptr<PlainBase> p;

}; // class RawFileReader



} // namespace ant

#endif // RAWFILEREADER_H
