#ifndef UNPACKER_H
#define UNPACKER_H

#include <TDataRecord.h>

#include <string>
#include <memory>
#include <list>



namespace ant {

class Unpacker {
public:

  class Module {
  public:
    virtual ~Module() = default;
    virtual bool OpenFile(const std::string& filename) = 0;
    virtual std::shared_ptr<TDataRecord> NextItem() = 0;
  };

  static std::unique_ptr<Module> Get(const std::string &filename);

  class Exception : public std::runtime_error {
    using std::runtime_error::runtime_error; // use base class constructor
  };

private:
  std::list< std::unique_ptr<Module> > modules;
  Unpacker();
};

} // namespace ant

#endif // UNPACKER_H
