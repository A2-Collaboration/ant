#ifndef UNPACKER_H
#define UNPACKER_H


#include <string>
#include <memory>
#include <list>

namespace ant {

class TDataRecord;

class Unpacker {

public:

  Unpacker() = delete;

  class Module {
  public:
    virtual ~Module() = default;
    virtual std::shared_ptr<TDataRecord> NextItem() = 0;
    virtual bool OpenFile(const std::string& filename) = 0;
  };

  // factory method to get a suitable unpacker module for the file
  static std::unique_ptr<Module> Get(const std::string &filename);

  class Exception : public std::runtime_error {
    using std::runtime_error::runtime_error; // use base class constructor
  };

};




} // namespace ant

#endif // UNPACKER_H
