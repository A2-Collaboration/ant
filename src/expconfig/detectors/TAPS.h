#ifndef DETECTORS_TAPS_H
#define DETECTORS_TAPS_H

#include "Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

#include "TVector2.h"
#include <limits>

namespace ant {
namespace expconfig {
namespace detector {


struct TAPS :
    ClusterDetector_t,
    UnpackerAcquConfig
{
  TAPS(
      bool cherenkovInstalled,
      bool useSensitiveChannels = false
      ) :
    ClusterDetector_t(Detector_t::Type_t::TAPS),
    CherenkovInstalled(cherenkovInstalled),
    UseSensitiveChannels(useSensitiveChannels)
  {
    // set up clusterelements pointers from derived class
    BuildClusterElements();
  }

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

  // ask derived class for actual elements, then this base class
  // can provide the ClusterDetector_t functions
  // (and makes some more checks)
  // we need references in order to modify the z-position for
  // the Cherenkov detector
  virtual std::vector<BaF2_Element_t>& GetBaF2Elements() const = 0;
  virtual std::vector<PbWO4_Element_t>& GetPbWO4Elements() const = 0;

private:
  void BuildClusterElements();
  // use another storage to make access to data performant
  std::vector<const ClusterDetector_t::Element_t*> clusterelements;
  bool CherenkovInstalled; // TAPS detectors moves downstream if Cherenkov installed
  bool UseSensitiveChannels; // Use sensitive channels as main integral
};


struct TAPS_2013 : TAPS {
  using TAPS::TAPS; // use constructor from base class

  virtual bool Matches(const THeaderInfo& headerInfo) const override;

protected:
  static std::vector<BaF2_Element_t>  BaF2_elements;
  static std::vector<PbWO4_Element_t> PbWO4_elements;
  virtual std::vector<BaF2_Element_t>& GetBaF2Elements() const override {
    return BaF2_elements;
  }
  virtual std::vector<PbWO4_Element_t>& GetPbWO4Elements() const override {
    return PbWO4_elements;
  }

}; // TAPS_2013

}}} // namespace ant::expconfig::detector

#endif // DETECTORS_TAPS_H
