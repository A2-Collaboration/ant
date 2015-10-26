#pragma once

#include "base/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

#include <cassert>

namespace ant {
namespace expconfig {
namespace detector {
struct EPT :
        TaggerDetector_t,
        UnpackerAcquConfig // EPT knows how to be filled from Acqu data
{

    virtual double GetPhotonEnergy(unsigned channel) const override {
        return BeamEnergy - elements[channel].ElectronEnergy;
    }
    virtual unsigned GetNChannels() const override {
        return elements.size();
    }
    virtual void SetIgnored(unsigned channel) override {
        elements[channel].Ignored = true;
    }
    virtual bool IsIgnored(unsigned channel) const override {
        return elements[channel].Ignored;
    }

    virtual bool TryGetChannelFromPhoton(double photonEnergy, unsigned& channel) const override;


    // for UnpackerAcquConfig
    virtual void BuildMappings(std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

protected:

    /// \todo have a look at ugcal?
    struct Element_t : TaggerDetector_t::Element_t {
        Element_t(
                unsigned channel,
                unsigned tdc,
                unsigned scaler,
                unsigned adc, // for Tagger, the ADC is least important
                double electronEnergy
                ) :
            TaggerDetector_t::Element_t(
                channel,
                electronEnergy
                ),
            TDC(tdc),
            Scaler(scaler),
            ADC(adc),
            Ignored(false)
        {}
        unsigned TDC;
        unsigned Scaler;
        unsigned ADC;
        bool Ignored;
    };

    EPT(double beamEnergy,
        const std::vector<Element_t>& elements_init) :
        TaggerDetector_t(
            Detector_t::Type_t::EPT,
            beamEnergy,
            3.2 // electronEnergyWidth
            ),
        elements(elements_init)
    {
        assert(elements.size()==47);
    }

    std::vector<Element_t> elements;
};

struct EPT_2014 : EPT {
    EPT_2014(double beamEnergy) :
        EPT(beamEnergy, elements_init)
    {}
    virtual bool Matches(const THeaderInfo& headerInfo) const override;
    static const std::vector<Element_t> elements_init;
};


}}} // namespace ant::expconfig::detector
