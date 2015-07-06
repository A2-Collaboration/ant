#ifndef BEAMTIMECONFIG_H
#define BEAMTIMECONFIG_H

#include <map>
#include <memory>

namespace ant {

class THeaderInfo;

class BeamtimeConfig
{
public:
  BeamtimeConfig() = delete;

  class Module {
  public:
    enum class Component_t {
      Trigger, Tagger, Moeller, CB, PID, MWPC, TAPS, TAPSVeto
    };

    enum class ChannelType_t {
      Time, Integral, IntegralShort,
      BitPattern, Scaler
    };

    struct LogicalChannel_t {
      Component_t Detector;
      ChannelType_t Type;
      unsigned LogicalChannel;
    };

    virtual ~Module() = default;

    virtual bool Matches(const THeaderInfo& header) = 0;

    using RawChannelMapping_t = std::map<unsigned, LogicalChannel_t>;
    virtual RawChannelMapping_t GetRawChannelMapping() = 0;
  };

  // factory method to get an unpacker
  static std::unique_ptr<Module> Get(const THeaderInfo& header);

  class Exception : public std::runtime_error {
    using std::runtime_error::runtime_error; // use base class constructor
  };
};

class BeamtimeConfigModule : public BeamtimeConfig::Module {};

} // namespace ant

#endif // BEAMTIMECONFIG_H
