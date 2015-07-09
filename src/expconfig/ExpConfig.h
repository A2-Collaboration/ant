#ifndef EXPCONFIG_H
#define EXPCONFIG_H

#include <map>
#include <memory>
#include <cstdint>


namespace ant {

class THeaderInfo;

enum class Detector_t : std::uint8_t {
  Trigger, Tagger, EPT, CB, PID, MWPC0, MWPC1,
  TAPS, TAPSVeto, Cherenkov, Moeller
};

enum class ChannelType_t : std::uint8_t {
  Timing, Integral, IntegralShort,
  BitPattern, Scaler, Counter
};

struct LogicalElement_t {
  Detector_t Detector;
  ChannelType_t Type;
  unsigned Element;
//  LogicalElement_t(
//      const Detector_t& detector,
//      const ChannelType_t& type,
//      unsigned element
//      ) :
//    Detector(detector),
//    Type(type),
//    Element(element)
//  {}
};

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
