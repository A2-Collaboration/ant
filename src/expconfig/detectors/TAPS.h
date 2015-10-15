#pragma once

#include "base/Detector_t.h"
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
    virtual unsigned GetNChannels() const override {
        return clusterelements.size();
    }
    virtual void SetIgnored(unsigned channel) override {
        clusterelements[channel]->Ignored = true;
    }
    virtual bool IsIgnored(unsigned channel) const override {
        return clusterelements[channel]->Ignored;
    }

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

    // for ClusterDetector_t
    virtual const ClusterDetector_t::Element_t* GetClusterElement(unsigned channel) const override {
        return clusterelements[channel];
    }

    /**
     * @brief GetRing returns the ring number inside the hexagonal structure
     * @param channel logical channel
     */
    unsigned GetRing(const unsigned channel) const;
    /**
     * @brief GetHexChannel returns the "usual" hexagonal channel as if four PbWO4s were one BaF2 hexagon
     * @param channel logical channel
     */
    unsigned GetHexChannel(const unsigned channel) const;

    /**
     * @brief IsPbWO4 returns true if given logical channel is PbWO4 element
     * @param channel logical channel
     * @return true if given logical channel is PbWO4 element
     * @note no boundary checks on channel performed
     */
    bool IsPbWO4(const unsigned channel) const;

    /**
     * @brief IsBaF2 returns true if given logical channel is BaF2 element
     * @param channel logical channel
     * @return true if given logical channel is BaF2 element
     * @note no boundary checks on channel performed
     */
    bool IsBaF2(const unsigned channel) const {
        return !IsPbWO4(channel);
    }

protected:

    // those values never changed over TAPS lifetime?!
    static constexpr unsigned NHexElements = 384;
    static constexpr unsigned NSectors = 6;
    static constexpr unsigned PbWO4PerHex = 4;

    // TAPS has BaF2 elements and PbWO4 elements

    struct TAPS_Element_t : ClusterDetector_t::Element_t {
        using ClusterDetector_t::Element_t::Element_t;
        bool Ignored = false;
    };

    struct BaF2_Element_t : TAPS_Element_t {
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
            TAPS_Element_t(
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

    struct PbWO4_Element_t : TAPS_Element_t {
        PbWO4_Element_t(
                unsigned channel,
                const TVector2& pos_xy,
                unsigned tdc,
                unsigned qdch,
                unsigned qdcl,
                const std::vector<unsigned>& neighbours // as last column since it's a messy list of numbers
                ) :
            TAPS_Element_t(
                channel,
                TVector3(pos_xy.X(), pos_xy.Y(), // z-component set by InitClusterElements()
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

    // given from derived class in constructor,
    // depending on the base class, the PbWO4_elements might be empty
    std::vector<BaF2_Element_t>  BaF2_elements;
    std::vector<PbWO4_Element_t> PbWO4_elements;

    // use another storage to make access to data performant
    void InitClusterElements();
    std::vector<TAPS_Element_t*> clusterelements;
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
