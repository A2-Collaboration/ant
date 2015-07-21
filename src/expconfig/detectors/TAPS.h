#ifndef DETECTORS_TAPS_H
#define DETECTORS_TAPS_H

#include "Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

#include "TVector2.h"
#include <limits>
#include <cassert>

#include <iostream>

namespace ant {
namespace expconfig {
namespace detector {


struct TAPS :
    ClusterDetector_t,
    UnpackerAcquConfig
{
  virtual TVector3 GetPosition(unsigned channel) const override {
    return clusterelements[channel]->Position;
  }

  virtual const ClusterDetector_t::Element_t* GetClusterElement(unsigned channel) const override {
    return clusterelements[channel];
  }

  // for UnpackerAcquConfig
  virtual void BuildMappings(
      std::vector<hit_mapping_t>&,
      std::vector<scaler_mapping_t>&) const override;

protected:

  // TAPS has BaF2 elements and PbWO4 elements

  struct BaF2_Element_t : ClusterDetector_t::Element_t {
    BaF2_Element_t(
        unsigned channel,
        const TVector2& pos_xy,
        unsigned tac,
        unsigned lg,
        unsigned sg,
        unsigned lgs,
        unsigned sgs,
        const std::vector<unsigned>& neighbours // as last column since it's a messy list of numbers
        ) :
      ClusterDetector_t::Element_t(
        channel,
        TVector3(pos_xy.X(), pos_xy.Y(), // z-component set by BuildClusterElements()
                 std::numeric_limits<double>::quiet_NaN()),
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
        const TVector2& pos_xy,
        unsigned tdc,
        unsigned qdch,
        unsigned qdcl,
        const std::vector<unsigned>& neighbours // as last column since it's a messy list of numbers
        ) :
      ClusterDetector_t::Element_t(
        channel,
        TVector3(pos_xy.X(), pos_xy.Y(), // z-component set by BuildClusterElements()
                 std::numeric_limits<double>::quiet_NaN()),
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

  TAPS(
      bool cherenkovInstalled,
      bool useSensitiveChannels,
      const std::vector<BaF2_Element_t>& BaF2s,
      const std::vector<PbWO4_Element_t>& PbWO4s
      ) :
    ClusterDetector_t(Detector_t::Type_t::TAPS),
    CherenkovInstalled(cherenkovInstalled),
    UseSensitiveChannels(useSensitiveChannels),
    BaF2_elements(BaF2s),
    PbWO4_elements(PbWO4s)
  {
    // init clusterelements from given BaF2/PbWO4 elements
    InitClusterElements();
  }


private:

  bool CherenkovInstalled; // TAPS detectors moves downstream if Cherenkov installed
  bool UseSensitiveChannels; // Use sensitive channels as main integral

  // given from derived class in constructor
  std::vector<BaF2_Element_t>  BaF2_elements;
  std::vector<PbWO4_Element_t> PbWO4_elements;

  // use another storage to make access to data performant
  void InitClusterElements();
  std::vector<const ClusterDetector_t::Element_t*> clusterelements;
};


struct TAPS_2013 : TAPS {
  TAPS_2013(
      bool cherenkovInstalled,
      bool useSensitiveChannels
      ) :
   TAPS(cherenkovInstalled, useSensitiveChannels,
        BaF2_elements_init, PbWO4_elements_init)
  {}

  virtual bool Matches(const THeaderInfo& headerInfo) const override;

private:
  const static std::vector<BaF2_Element_t>  BaF2_elements_init;
  const static std::vector<PbWO4_Element_t> PbWO4_elements_init;

}; // TAPS_2013

}}} // namespace ant::expconfig::detector

#endif // DETECTORS_TAPS_H
