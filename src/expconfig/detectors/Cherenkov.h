#pragma once

#include "base/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {
struct Cherenkov :
        Detector_t,
        UnpackerAcquConfig // Cherenkov knows how to be filled from Acqu data
{
    Cherenkov() :
        Detector_t(Type_t::Cherenkov),
        /// \todo think of a better position, fix cherenkov ADC input channel number (17 is wrong)
        element(0, {std_ext::NaN,std_ext::NaN,std_ext::NaN}, 17)
    {}

    virtual vec3 GetPosition(unsigned) const override {
        return element.Position;
    }
    virtual unsigned GetNChannels() const override {
        return 1;
    }
    virtual void SetElementFlags(unsigned, const ElementFlags_t& flags) override {
        element.Flags |= flags;
    }
    virtual const ElementFlags_t& GetElementFlags(unsigned) const override {
        return element.Flags;
    }

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>& hit_mappings,
            std::vector<scaler_mapping_t>&) const override
    {
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  element.ADC);
    }

protected:

    struct Element_t : Detector_t::Element_t {
         Element_t(unsigned channel, const vec3& position, unsigned adc) :
             Detector_t::Element_t(channel, position), ADC(adc) {}
         unsigned ADC;
    };

    Element_t element;

};

}}} // namespace ant::expconfig::detector
