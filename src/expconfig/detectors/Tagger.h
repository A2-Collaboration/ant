#pragma once

#include "base/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

#include <cassert>

namespace ant {
namespace expconfig {
namespace detector {
struct Tagger :
        TaggerDetector_t,
        UnpackerAcquConfig // EPT knows how to be filled from Acqu data
{

    virtual double GetPhotonEnergy(unsigned channel) const override {
        return BeamEnergy - elements[channel].ElectronEnergy;
    }
    virtual ant::TaggerDetector_t::taggeff_t GetTaggEff(unsigned channel) const override
    {
        return elements.at(channel).TaggEff;
    }
    virtual void SetTaggEff(unsigned channel, const ant::TaggerDetector_t::taggeff_t& taggeff) override
    {
        elements.at(channel).TaggEff = taggeff;
    }
    virtual double GetPhotonEnergyWidth(unsigned) const override {
        return 3.2; /// \todo Find better approximation for tagger here!
    }
    virtual unsigned GetNChannels() const override {
        return unsigned(elements.size());
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

    Tagger(double beamEnergy,
        const std::vector<Element_t>& elements_init) :
        TaggerDetector_t(
            Detector_t::Type_t::Tagger,
            beamEnergy
            ),
        elements(elements_init)
    {
        assert(elements.size()==352);
    }

    std::vector<Element_t> elements;
};

struct Tagger_2007 : Tagger {
    Tagger_2007() :
        Tagger(1508.0, elements_init)
    {}
    static const std::vector<Element_t> elements_init;
};

struct Tagger_2010_03 : Tagger {
    Tagger_2010_03() :
        Tagger(450.0, elements_init)
    {}
    static const std::vector<Element_t> elements_init;
};

struct Tagger_2015 : Tagger {
    Tagger_2015() :
        Tagger(450.0, elements_init)
    {}
    static const std::vector<Element_t> elements_init;
};


}}} // namespace ant::expconfig::detector
