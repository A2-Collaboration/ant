#ifndef DETECTORS_TAPS_H
#define DETECTORS_TAPS_H

#include "Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {


struct TAPS :
    ClusterDetector_t,
    UnpackerAcquConfig
{
  TAPS(bool cherenkovInstalled) :
    ClusterDetector_t(Detector_t::Type_t::TAPS),
    CherenkovInstalled(cherenkovInstalled)
  {}

protected:
  // TAPS has BaF2 elements and PbWO4 elements

  struct BaF2_Element_t : ClusterDetector_t::Element_t {
    BaF2_Element_t(
        unsigned channel,
        const TVector3& position,
        unsigned tac,
        unsigned lg,
        unsigned sg,
        unsigned lgs,
        unsigned sgs,
        const std::vector<unsigned>& neighbours // as last column since it's a messy list of numbers
        ) :
      ClusterDetector_t::Element_t(
        channel,
        position,
        neighbours,
        3.4 /// \todo use best value from S. Lohse diploma thesis?
        ),
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

  struct PbWO4_Element_t : ClusterDetector_t::Element_t {
    PbWO4_Element_t(
        unsigned channel,
        const TVector3& position,
        unsigned tdc,
        unsigned qdch,
        unsigned qdcl,
        const std::vector<unsigned>& neighbours // as last column since it's a messy list of numbers
        ) :
      ClusterDetector_t::Element_t(
        channel,
        position,
        neighbours,
        2.2 /// \todo use best value from S. Lohse diploma thesis?
        ),
      TDC(tdc),
      QDCH(qdch),
      QDCL(qdcl)
    {}
    unsigned TDC;  // timing
    unsigned QDCH; // integral
    unsigned QDCL; // integral, sensitive
  };

  bool CherenkovInstalled; // TAPS detectors moves downstream if Cherenkov installed
};


struct TAPS_2013 : TAPS {
  using TAPS::TAPS; // use constructor from base class


  virtual bool Matches(const THeaderInfo& headerInfo) const override;

  virtual TVector3 GetPosition(unsigned channel) const override {
    /// \todo Implement proper position stuff for TAPS
    return TVector3();
  }

  virtual const ClusterDetector_t::Element_t* GetClusterElement(unsigned channel) const override {
    return nullptr;
  }

  // for UnpackerAcquConfig
  virtual void BuildMappings(
      std::vector<hit_mapping_t>&,
      std::vector<scaler_mapping_t>&) const override;

protected:
  static const std::vector<BaF2_Element_t>  BaF2_elements;
  static const std::vector<PbWO4_Element_t> PbWO4_elements;


}; // TAPS_2013

}}} // namespace ant::expconfig::detector

#endif // DETECTORS_TAPS_H
