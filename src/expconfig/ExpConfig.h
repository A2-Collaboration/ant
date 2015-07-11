#ifndef EXPCONFIG_H
#define EXPCONFIG_H

#include <map>
#include <memory>
#include <cstdint>


namespace ant {

class THeaderInfo;

struct Channel_t {
  enum class Type_t : std::uint8_t {
    Timing,
    Integral, IntegralShort,
    IntegralAlternate, IntegralShortAlternate,
    BitPattern, Scaler, Counter
  };
  static bool IsIntegral(const Type_t t);
};

struct Detector_t {
  enum class Type_t : std::uint8_t {
    Trigger, Tagger, EPT, CB, PID, MWPC0, MWPC1,
    TAPS, TAPSVeto, Cherenkov, Moeller
  };
  const Type_t Type;

  // Element_t is the minimum information,
  // derived classes may extend this class
  struct Element_t {
    struct Position_t {
      double X, Y, Z;
    };
    Element_t(unsigned channel, const Position_t& position) :
      Channel(channel),
      Position(position)
    {}
    unsigned Channel; // unique within Detector for all time!
    Position_t Position;
  };

  virtual ~Detector_t() = default;

protected:
  Detector_t(const Type_t& type) :
    Type(type) {}
  Detector_t(const Detector_t&) = delete; // disable copy
};



struct LogicalChannel_t {
  Detector_t::Type_t Detector;
  Channel_t::Type_t Type;
  unsigned Channel;
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
