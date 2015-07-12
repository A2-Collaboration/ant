#ifndef EXPCONFIG_DETECTOR_T_H
#define EXPCONFIG_DETECTOR_T_H

#include <string>
#include <vector>

namespace ant {

/**
 * @brief The Detector_t struct is the minimal base class for all detectors
 */
struct Detector_t {
  enum class Type_t : std::uint8_t {
    Trigger, Tagger, EPT, CB, PID, MWPC0, MWPC1,
    TAPS, TAPSVeto, Cherenkov, Moeller
  };
  const Type_t Type;
  static std::string ToString(const Type_t& type);

  // Element_t is the minimum information,
  // derived classes may extend this
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

/**
 * @brief The Channel_t struct enumerates the different channel types
 */
struct Channel_t {
  enum class Type_t : std::uint8_t {
    Timing,
    Integral, IntegralShort,
    IntegralAlternate, IntegralShortAlternate,
    BitPattern, Scaler, Counter
  };
  static bool IsIntegral(const Type_t& t);
  static std::string ToString(const Type_t& type);
};

/**
 * @brief The LogicalChannel_t struct uniquely identifies detector element within setup
 */
struct LogicalChannel_t {
  Detector_t::Type_t Detector;
  Channel_t::Type_t Type;
  unsigned Channel;
};

struct ClusterDetector_t : Detector_t {

  virtual std::vector<unsigned> GetNeighbours(unsigned channel) const = 0;
  virtual double GetMoliereRadius(unsigned channel) const = 0;


protected:
  ClusterDetector_t(const Type_t& type) :
    Detector_t(type) {}
};

} // namespace ant

#endif // EXPCONFIG_DETECTOR_T_H
