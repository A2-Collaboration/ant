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

    virtual vec3 GetPosition(unsigned channel) const override {
        return elements.at(channel)->Position;
    }
    virtual unsigned GetNChannels() const override {
        return elements.size();
    }
    virtual void SetElementFlags(unsigned channel, const ElementFlags_t& flags) override {
        elements[channel]->Flags |= flags;
    }
    virtual const ElementFlags_t& GetElementFlags(unsigned channel) const override {
        return elements[channel]->Flags;
    }

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

    /**
     * @brief GetZPosition distance of front face from center of target
     * @return distance in centimeters
     */
    virtual double GetZPosition() const;

    /**
     * @brief Get a rough diameter containing the veto element
     * @return diameter in cm
     * @note Roughly estimated by looking at a dummy TAPS element
     */
    double GetElementDiameter() const { return 7.0; }

    bool IsPbWO4(const unsigned channel) const;

protected:

    static constexpr unsigned NSectors = 6;

    using TAPSVeto_Element_t = Detector_t::Element_t;

    struct BaF2_Element_t : TAPSVeto_Element_t {
        BaF2_Element_t(
                unsigned channel,
                const vec2& pos_xy,
                unsigned tac,
                unsigned lgs // only sensitive for BaF2 Veto
                ) :
            TAPSVeto_Element_t(
                channel,
                vec3(pos_xy.x, pos_xy.y, // z-component set by InitElements()
                     std::numeric_limits<double>::quiet_NaN())
                ),
            TAC(tac),
            LGS(lgs)
        {}
        unsigned TAC;
        unsigned LGS;
    };

    struct PbWO4_Element_t : TAPSVeto_Element_t {
      PbWO4_Element_t(
          unsigned channel,
          const vec2& pos_xy,
          unsigned tdc,
          unsigned qdch,
          unsigned qdcl
          ) :
        TAPSVeto_Element_t(
          channel,
          vec3(pos_xy.x, pos_xy.y, // z-component set by InitClusterElements()
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
            bool pizzaInstalled,
            const std::vector<BaF2_Element_t>& BaF2s,
            const std::vector<PbWO4_Element_t>& PbWO4s) :
        Detector_t(Detector_t::Type_t::TAPSVeto),
        CherenkovInstalled(cherenkovInstalled),
        PizzaInstalled(pizzaInstalled),
        BaF2_elements(BaF2s),
        PbWO4_elements(PbWO4s)
    {
        // init clusterelements from given BaF2/PbWO4 elements
        InitElements();
    }



private:

    bool CherenkovInstalled;  // TAPS detectors moves downstream if Cherenkov installed
    bool PizzaInstalled;  // TAPS moves downstream as well if the Pizza detector is installed

    // given from derived class in constructor
    std::vector<BaF2_Element_t>  BaF2_elements;
    std::vector<PbWO4_Element_t> PbWO4_elements;

    void InitElements();
    std::vector<TAPSVeto_Element_t*> elements;

};


struct TAPSVeto_2014 : TAPSVeto {
    TAPSVeto_2014(
            bool cherenkovInstalled,
            bool pizzaInstalled  // Pizza is only available starting with this configuration
            ) :
        TAPSVeto(cherenkovInstalled, pizzaInstalled,
                 BaF2_elements_init, PbWO4_elements_init)
    {}

protected:
    TAPSVeto_2014(
            bool cherenkovInstalled,
            bool pizzaInstalled,
            const std::vector<BaF2_Element_t>& BaF2s
            ) :
        TAPSVeto(cherenkovInstalled, pizzaInstalled,
                 BaF2s, PbWO4_elements_init)
    {}

private:
    const static std::vector<BaF2_Element_t>  BaF2_elements_init;
    const static std::vector<PbWO4_Element_t> PbWO4_elements_init;
};

struct TAPSVeto_2013_11 : TAPSVeto {
    TAPSVeto_2013_11(
            bool cherenkovInstalled
            ) :
        TAPSVeto(cherenkovInstalled, false,
                 BaF2_elements_init, PbWO4_elements_init)
    {}

protected:
    TAPSVeto_2013_11(
            bool cherenkovInstalled,
            const std::vector<BaF2_Element_t>& BaF2s
            ) :
        TAPSVeto(cherenkovInstalled, false,
                 BaF2s, PbWO4_elements_init)
    {}

private:
    const static std::vector<BaF2_Element_t>  BaF2_elements_init;
    const static std::vector<PbWO4_Element_t> PbWO4_elements_init;
};

struct TAPSVeto_2009_03 : TAPSVeto {
    TAPSVeto_2009_03(
            bool cherenkovInstalled
            ) :
        TAPSVeto(cherenkovInstalled, false,
                 BaF2_elements_init, PbWO4_elements_init)
    {}

protected:
    TAPSVeto_2009_03(
            bool cherenkovInstalled,
            const std::vector<BaF2_Element_t>& BaF2s
            ) :
        TAPSVeto(cherenkovInstalled, false,
                 BaF2s, PbWO4_elements_init)
    {}

private:
    const static std::vector<BaF2_Element_t>  BaF2_elements_init;
    const static std::vector<PbWO4_Element_t> PbWO4_elements_init;
};

struct TAPSVeto_2007: TAPSVeto {
    TAPSVeto_2007(
            bool cherenkovInstalled
            ) :
        TAPSVeto(cherenkovInstalled, false,
                 BaF2_elements_init, {})
    {}

protected:
    TAPSVeto_2007(
            bool cherenkovInstalled,
            const std::vector<BaF2_Element_t>& BaF2s
            ) :
        TAPSVeto(cherenkovInstalled, false,
                 BaF2s, {})
    {}

private:
    const static std::vector<BaF2_Element_t>  BaF2_elements_init;
};


}}} // namespace ant::expconfig::detector
