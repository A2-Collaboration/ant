#include "Detector_t.h"

#include <stdexcept>

using namespace std;
using namespace ant;

bool Channel_t::IsIntegral(const Channel_t::Type_t& t) {
  switch(t) {
  case Type_t::Integral:
  case Type_t::IntegralShort:
  case Type_t::IntegralAlternate:
  case Type_t::IntegralShortAlternate:
    return true;
  default:
    return false;
  }
}

const char *Channel_t::ToString(const Type_t& type)
{
  switch(type) {
  case Channel_t::Type_t::BitPattern:
    return "BitPattern";
  case Channel_t::Type_t::Counter:
    return "Counter";
  case Channel_t::Type_t::Integral:
    return "Integral";
  case Channel_t::Type_t::IntegralAlternate:
    return "IntegralAlternate";
  case Channel_t::Type_t::IntegralShort:
    return "IntegralShort";
  case Channel_t::Type_t::IntegralShortAlternate:
    return "IntegralShortAlternate";
  case Channel_t::Type_t::Scaler:
    return "Scaler";
  case Channel_t::Type_t::Timing:
    return "Timing";
  }
  throw runtime_error("Not implemented");
}

const char* Detector_t::ToString(const Type_t &type)
{
  switch(type) {
  case Detector_t::Type_t::CB :
    return "CB";
  case Detector_t::Type_t::Cherenkov:
    return "Cherenkov";
  case Detector_t::Type_t::EPT:
    return "EPT";
  case Detector_t::Type_t::Moeller:
    return "Moeller";
  case Detector_t::Type_t::MWPC0:
    return "MWPC0";
  case Detector_t::Type_t::MWPC1:
    return "MWPC1";
  case Detector_t::Type_t::PID:
    return "PID";
  case Detector_t::Type_t::Tagger:
    return "Tagger";
  case Detector_t::Type_t::TAPS:
    return "TAPS";
  case Detector_t::Type_t::TAPSVeto:
    return "TAPSVeto";
  case Detector_t::Type_t::Trigger:
    return "Trigger";
  }
  throw runtime_error("Not implemented");
}
