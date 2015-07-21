#ifndef DETECTORS_EPT_H
#define DETECTORS_EPT_H

#include "Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {
struct EPT :
        TaggerDetector_t,
        UnpackerAcquConfig // EPT knows how to be filled from Acqu data
{

    EPT(double beamEnergy) :
        TaggerDetector_t(Detector_t::Type_t::EPT, beamEnergy) {}

    virtual bool Matches(const THeaderInfo&) const override {
        // always match, since EPT never changed over A2's lifetime
        return true;
    }

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

protected:
    struct Element_t : TaggerDetector_t::Element_t {
        Element_t(
                unsigned channel,
                double electronEnergy,
                unsigned adc,
                unsigned tdc,
                unsigned scaler
                ) :
            TaggerDetector_t::Element_t(
                channel,
                electronEnergy
                ),
            ADC(adc),
            TDC(tdc),
            Scaler(scaler)
        {}
        unsigned ADC;
        unsigned TDC;
        unsigned Scaler;
    };
    static const std::vector<Element_t> elements;
};


}}} // namespace ant::expconfig::detector



#endif // DETECTORS_EPT_H
