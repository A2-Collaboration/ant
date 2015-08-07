#pragma once

#include "base/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {
struct TAPSVeto :
        Detector_t,
        UnpackerAcquConfig // TAPSVeto knows how to be filled from Acqu data
{

    virtual TVector3 GetPosition(unsigned channel) const override {
        return elements.at(channel)->Position;
    }
    virtual unsigned GetNChannels() const override {
        return elements.size();
    }

    virtual bool Matches(const THeaderInfo&) const override {
        // always match, since TAPSVeto never changed over A2's lifetime
        return true;
    }

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

    /**
     * @brief Get a radius containing the veto element. Used for candidate building
     * @return radius in cm
     * @note Roughly estimated by looking at a dummy TAPS element
     */
    double GetElementRadius() const { return 7.0; }

protected:

    struct BaF2_Element_t : Detector_t::Element_t {
        BaF2_Element_t(
                unsigned channel,
                const TVector2& pos_xy,
                unsigned tac,
                unsigned lgs // only sensitive for BaF2 Veto
                ) :
            Detector_t::Element_t(
                channel,
                TVector3(pos_xy.X(), pos_xy.Y(), // z-component set by InitElements()
                         std::numeric_limits<double>::quiet_NaN())
                ),
            TAC(tac),
            LGS(lgs)
        {}
        unsigned TAC;
        unsigned LGS;
    };

    struct PbWO4_Element_t : Detector_t::Element_t {
      PbWO4_Element_t(
          unsigned channel,
          const TVector2& pos_xy,
          unsigned tdc,
          unsigned qdch,
          unsigned qdcl
          ) :
        Detector_t::Element_t(
          channel,
          TVector3(pos_xy.X(), pos_xy.Y(), // z-component set by InitClusterElements()
                   std::numeric_limits<double>::quiet_NaN())
          ),
        TDC(tdc),
        QDCH(qdch),
        QDCL(qdcl)
      {}
      unsigned TDC;  // timing
      unsigned QDCH; // integral
      unsigned QDCL; // integral, sensitive?
    };

    TAPSVeto(
            bool cherenkovInstalled,
            const std::vector<BaF2_Element_t>& BaF2s,
            const std::vector<PbWO4_Element_t>& PbWO4s) :
        Detector_t(Detector_t::Type_t::TAPSVeto),
        CherenkovInstalled(cherenkovInstalled),
        BaF2_elements(BaF2s),
        PbWO4_elements(PbWO4s)
    {
        // init clusterelements from given BaF2/PbWO4 elements
        InitElements();
    }



private:

    bool CherenkovInstalled; // TAPS detectors moves downstream if Cherenkov installed

    // given from derived class in constructor
    std::vector<BaF2_Element_t>  BaF2_elements;
    std::vector<PbWO4_Element_t> PbWO4_elements;

    void InitElements();
    std::vector<const Detector_t::Element_t*> elements;

};


struct TAPSVeto_2013 : TAPSVeto {
    TAPSVeto_2013(
            bool cherenkovInstalled
            ) :
        TAPSVeto(cherenkovInstalled,
                 BaF2_elements_init, PbWO4_elements_init)
    {}

    virtual bool Matches(const THeaderInfo& headerInfo) const override;

protected:
    TAPSVeto_2013(
            bool cherenkovInstalled,
            const std::vector<BaF2_Element_t>& BaF2s
            ) :
        TAPSVeto(cherenkovInstalled,
                 BaF2s, PbWO4_elements_init)
    {}

private:
    const static std::vector<BaF2_Element_t>  BaF2_elements_init;
    const static std::vector<PbWO4_Element_t> PbWO4_elements_init;
};

struct TAPSVeto_2014 : TAPSVeto_2013 {
    TAPSVeto_2014(
            bool cherenkovInstalled
            ) :
        TAPSVeto_2013(cherenkovInstalled, BaF2_elements_init)
    {}

    virtual bool Matches(const THeaderInfo& headerInfo) const override;

private:
    const static std::vector<BaF2_Element_t>  BaF2_elements_init;
}; // TAPSVeto_2014

}}} // namespace ant::expconfig::detector
