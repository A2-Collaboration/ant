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
    return positions.at(channel);
  }

  virtual bool Matches(const THeaderInfo&) const override {
    // always match, since CB never changed over A2's lifetime
    return true;
  }

  // for UnpackerAcquConfig
  virtual void BuildMappings(
      std::vector<hit_mapping_t>&,
      std::vector<scaler_mapping_t>&) const override;

  virtual double GetMoliereRadius(unsigned) const override {
    return 4.5;
  }

  virtual std::vector<unsigned> GetNeighbours(unsigned) const override {
    return {};
  }

protected:
  struct CBElement_t : Element_t {
    CBElement_t(
        unsigned channel,
        const TVector3& position,
        unsigned adc,
        unsigned tdc
        ) :
      Element_t(channel, position), // init fields
      ADC(adc),
      TDC(tdc)
    {}
    unsigned ADC;
    unsigned TDC;
  };
  static const std::vector<CBElement_t> elements;
  static const std::map<unsigned, TVector3> positions;
  static std::vector<CBElement_t> initElements();
  static std::map<unsigned, TVector3> initPositions();
};


}}} // namespace ant::expconfig::detector



#endif // DETECTORS_CB_H
