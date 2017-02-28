#pragma once

#include "base/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {
struct CB :
        ClusterDetector_t,
        UnpackerAcquConfig // CB knows how to be filled from Acqu data
{

    CB();

    virtual vec3 GetPosition(unsigned channel) const override {
        return elements[channel].Position;
    }

    virtual double GetInnerRadius() const {
        return 25.4; // the famous 10 inch
    }

    virtual unsigned GetNChannels() const override {
        return elements.size();
    }

    virtual void SetElementFlags(unsigned channel, const ElementFlags_t& flags) override;
    virtual const ElementFlags_t& GetElementFlags(unsigned channel) const override {
        return elements.at(channel).Flags;
    }
    virtual bool IsHole(unsigned channel) const;

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

    // for ClusterDetector_t
    virtual const ClusterDetector_t::Element_t* GetClusterElement(unsigned channel) const override {
        return std::addressof(elements[channel]);
    }

protected:
    struct Element_t : ClusterDetector_t::Element_t {
        Element_t(
                unsigned channel,
                const vec3& position,
                unsigned adc,
                unsigned tdc,
                const std::vector<unsigned>& neighbours
                ) :
            ClusterDetector_t::Element_t(
                channel,
                position,
                neighbours,
                4.8, /// \todo use best value from S. Lohse diploma thesis?
                13.3, // critical energy
                2.588 // radiation length
                ),
            ADC(adc),
            TDC(tdc)
        {}
        unsigned ADC;
        unsigned TDC;
    };
    static const std::vector<Element_t> elements_init;
    std::vector<Element_t> elements;

    void SetTouchesHoleOfNeighbours(unsigned hole);
};


}}} // namespace ant::expconfig::detector
