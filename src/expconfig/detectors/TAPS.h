#ifndef DETECTORS_TAPS_H
#define DETECTORS_TAPS_H

#include "ExpConfig.h"


namespace ant {
namespace expconfig {
namespace detector {


struct TAPS : Detector_t {
  TAPS() : Detector_t(Detector_t::Type_t::TAPS) {}
};


}}} // namespace ant::expconfig::detector

#endif // DETECTORS_TAPS_H
