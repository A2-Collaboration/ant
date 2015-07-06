#ifndef UNPACKER_H
#define UNPACKER_H

#include <TDataRecord.h>

#include <string>
#include <memory>
#include <list>

namespace ant {

class THeaderInfo;

class Unpacker {
public:

  Unpacker() = delete;

  class Module {
  public:
    virtual ~Module() = default;
    virtual std::shared_ptr<TDataRecord> NextItem() = 0;
    virtual bool OpenFile(const std::string& filename) = 0;
  };

  // factory method to get an unpacker module
  static std::unique_ptr<Module> Get(const std::string &filename);

  template<class Module>
  class Config {
  public:
    virtual bool Matches(const THeaderInfo& headerInfo) = 0;
    static std::unique_ptr< Unpacker::Config<Module> > Get(const THeaderInfo& headerInfo);
  };

  class Exception : public std::runtime_error {
    using std::runtime_error::runtime_error; // use base class constructor
  };

};




} // namespace ant

#endif // UNPACKER_H
