#ifndef EXPCONFIG_H
#define EXPCONFIG_H

#include <memory>
#include "base/printable.h"
#include <list>
#include "Detector_t.h"

namespace ant {

class THeaderInfo;
class CalibrationApply_traits;
class CalibrationUpdate_traits;
class Detector_t;

class ExpConfig
{
public:
  ExpConfig() = delete; // this class is more a wrapper for handling the config

  // all configs have a common base and match via THeaderInfo
  class Base {
  public:
    virtual ~Base() = default;
    virtual bool Matches(const THeaderInfo& header) const = 0;
  };

  // the ExpConfig::Module provides general information about the experiment
  class Module : public virtual Base {
  public:
    static std::unique_ptr<Module> Get(const THeaderInfo& header);
  };

  // in order to run the Reconstruction,
  // the following methods are needed
  class Reconstruct : public virtual Base {
  public:
    virtual std::list< std::shared_ptr< CalibrationApply_traits > > GetCalibrations() const = 0;
    virtual std::list< std::shared_ptr< CalibrationUpdate_traits > > GetUpdateables() const = 0;
    static std::unique_ptr<Reconstruct> Get(const THeaderInfo& header);
  };

  // each unpacker has its own config,
  // we enforce this by templates
  template<class T>
  class Unpacker : public virtual Base {
  public:
    static std::unique_ptr<T> Get(const THeaderInfo& header);
  };

  class Exception : public std::runtime_error {
    using std::runtime_error::runtime_error; // use base class constructor
  };
};

} // namespace ant



#endif // EXPCONFIG_H
