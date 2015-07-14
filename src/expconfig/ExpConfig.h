#ifndef EXPCONFIG_H
#define EXPCONFIG_H

#include <map>
#include <memory>
#include <cstdint>
#include "base/printable.h"

#include "Detector_t.h"

namespace ant {

class THeaderInfo;
class CalibrationApply_traits;

class ExpConfig
{
public:
  ExpConfig() = delete;

  // all configs have a common base and match via THeaderInfo
  class Base {
  public:
    virtual ~Base() = default;
    virtual bool Matches(const THeaderInfo& header) const = 0;
  };

  // the ExpConfig::Module defines methods for unpacker-independent configs
  class Module : public virtual Base {

    virtual std::vector< std::unique_ptr< CalibrationApply_traits > > GetCalibrations() const = 0;
  };

  // each unpacker has its own config,
  // we enforce this by templates
  template<class T>
  class Unpacker : public virtual Base {
  public:
    static std::unique_ptr<T> Get(const THeaderInfo& header);
  };

  // factory method to get a config for already unpacked data
  static std::unique_ptr<Module> Get(const THeaderInfo& header);

  class Exception : public std::runtime_error {
    using std::runtime_error::runtime_error; // use base class constructor
  };
};

} // namespace ant



#endif // EXPCONFIG_H
