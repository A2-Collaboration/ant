#ifndef DETECTORS_CB_H
#define DETECTORS_CB_H

#include "ExpConfig.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {
struct CB : Detector_t, UnpackerAcquConfig {

  CB() : Detector_t(Detector_t::Type_t::CB) {}



  virtual bool Matches(const THeaderInfo &header) const override {
    return true;
  }

  virtual void BuildMappings(std::vector<hit_mapping_t> &hit_mappings,
                             std::vector<scaler_mapping_t> &) {
    // CB has only hit_mappings to add, no scalers

  }

private:
  struct CBElement_t : Element_t {
    CBElement_t(unsigned channel,
                const Position_t& position,
                unsigned adc,
                unsigned tdc) :
      Element_t(channel, position), // init fields
      ADC(adc),
      TDC(tdc)
    {}
    unsigned ADC;
    unsigned TDC;
  };
  static const std::vector<CBElement_t> elements;



};


}}} // namespace ant::expconfig::detector



#endif // DETECTORS_CB_H
