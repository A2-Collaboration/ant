#ifndef RAWFILEREADER_H
#define RAWFILEREADER_H

#include <fstream>
#include <string>
#include <memory>
#include <cstdint>

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
  
  explicit RawFileReader(const std::string& filename) 
  {
    if(std_ext::string_ends_with(filename, ".xz")) {
      p = std::unique_ptr<Plain>(new XZ(filename));
    }
    else {
      p = std::unique_ptr<Plain>(new Plain(filename));
    }
  }
  
  // methods/operators for testing the state
  //bool is_open() const { return file.is_open(); }
  //bool operator!() const { return file.operator!(); }
  explicit operator bool() const { 
    return p->operator bool();
  }
  
  // read
  void read(char* s, std::streamsize n) { 
    p->read(s,n); 
  }
  
  std::streamsize gcount() const {
    return p->gcount();
  }
  
  class Exception : public std::runtime_error {
    using std::runtime_error::runtime_error; // use base class constructor
  };
  
private:
  
  /**
   * @brief The Plain class
   * 
   * Plain and simple filereader, encapsulating an ifstream 
   * and providing a base class for the more complicated decompression classes
   * 
   * Only the really needed methods are exported
   */
  class Plain {
  public:
    explicit Plain(const std::string& filename) 
      : file(filename.c_str(), std::ios::binary)
    {}
    
    virtual ~Plain() {} // derived classes may want to 
    
    // methods/operators for testing the state
    //bool is_open() const { return file.is_open(); }
    //bool operator!() const { return file.operator!(); }
    virtual explicit operator bool() const { 
      return file.operator bool(); 
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
  class XZ : public Plain {
  public:  
    
    XZ(const std::string& filename);
    
    virtual ~XZ();
    
    virtual explicit operator bool() const { 
      return Plain::operator bool() && !decompressFailed; 
    }
    
    virtual void read(char *s, std::streamsize n);
    
    virtual std::streamsize gcount() const {
      return gcount_;
    }
   
  private:
    bool decompressFailed;
    std::streamsize gcount_;
    
#ifndef RAWFILEREADER_H_IMPL
    struct lzma_stream;
#endif   
    lzma_stream* strm;    
    void init_decoder();
    
    void cleanup();
    
    
    
  }; // class RawFileReader::XZ
  
  
  // private stuff for RawFileReader
  std::unique_ptr<Plain> p;

  
}; // class RawFileReader



} // namespace ant

#endif // RAWFILEREADER_H
