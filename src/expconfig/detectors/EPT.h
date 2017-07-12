#pragma once

#include "base/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

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
    virtual void SetElementFlags(unsigned channel, const ElementFlags_t& flags) override {
        elements[channel].Flags |= flags;
    }
    virtual const ElementFlags_t& GetElementFlags(unsigned channel) const override {
        return elements[channel].Flags;
    }

    // for UnpackerAcquConfig
    virtual void BuildMappings(std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

    const static std::string ScalerName;

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
            ADC(adc)
        {}
        unsigned TDC;
        unsigned Scaler;
        unsigned ADC;
    };

    EPT(double beamEnergy,
        const std::vector<Element_t>& elements_init);

    std::vector<Element_t> elements;
};

struct EPT_2014 : EPT {
    EPT_2014(double beamEnergy) :
        EPT(beamEnergy, elements_init)
    {}
    static const std::vector<Element_t> elements_init;
};


}}} // namespace ant::expconfig::detector
