#ifndef DETECTORS_TAPS_H
#define DETECTORS_TAPS_H

#include "ExpConfig.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {


struct TAPS :
    Detector_t,
    UnpackerAcquConfig
{
  TAPS(bool cherenkovInstalled) :
    Detector_t(Detector_t::Type_t::TAPS),
    CherenkovInstalled(cherenkovInstalled)
  {}

protected:
  struct TAPSElement_t : Element_t {
    TAPSElement_t(unsigned channel,
                const Position_t& position,
                unsigned tac,
                unsigned lg,
                unsigned sg,
                unsigned lgs,
                unsigned sgs
                ) :
      Element_t(channel, position), // init fields
      TAC(tac),
      LG(lg),
      SG(sg),
      LGS(lgs),
      SGS(sgs)
    {}
    unsigned TAC; // timing
    unsigned LG;  // integral, long gate
    unsigned SG;  // integral, short gate
    unsigned LGS; // integral, long  gate, high gain
    unsigned SGS; // integral, short gate, high gain

  };

  bool CherenkovInstalled; // TAPS detectors moves downstream if Cherenkov installed

};


struct TAPS_2013 : TAPS {
  using TAPS::TAPS; // use constructor from base class


  virtual bool Matches(const THeaderInfo& headerInfo) const override;

  // for UnpackerAcquConfig
  virtual void BuildMappings(
      std::vector<hit_mapping_t>&,
      std::vector<scaler_mapping_t>&) const override;
};

}}} // namespace ant::expconfig::detector

#endif // DETECTORS_TAPS_H
