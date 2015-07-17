#ifndef DETECTORS_CB_H
#define DETECTORS_CB_H

#include "Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {
struct CB :
    ClusterDetector_t,
    UnpackerAcquConfig // CB knows how to be filled from Acqu data
{

  CB() : ClusterDetector_t(Detector_t::Type_t::CB) {}

  virtual TVector3 GetPosition(unsigned channel) const override {
    return elements[channel].Position;
  }

  virtual bool Matches(const THeaderInfo&) const override {
    // always match, since CB never changed over A2's lifetime
    return true;
  }

  // for UnpackerAcquConfig
  virtual void BuildMappings(
      std::vector<hit_mapping_t>&,
      std::vector<scaler_mapping_t>&) const override;

  virtual const ClusterElement_t* GetClusterElement(unsigned channel) const override {
    return std::addressof(elements[channel]);
  }

protected:
  struct CBElement_t : ClusterElement_t {
    CBElement_t(
        unsigned channel,
        const TVector3& position,
        unsigned adc,
        unsigned tdc,
        const std::vector<unsigned>& neighbours
        ) :
      ClusterElement_t(channel, position, neighbours, 4.0), // all NaI elements have 4.0 as MoliereRadius
      ADC(adc),
      TDC(tdc)
    {}
    unsigned ADC;
    unsigned TDC;
  };
  static const std::vector<CBElement_t> elements;
  static const std::vector<unsigned> holes;

  static std::vector<unsigned> initHoles();
};


}}} // namespace ant::expconfig::detector



#endif // DETECTORS_CB_H
