#pragma once

#include "expconfig/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {
struct PID :
        Detector_t,
        UnpackerAcquConfig // PID knows how to be filled from Acqu data
{


    virtual TVector3 GetPosition(unsigned channel) const override {
        return elements[channel].Position;
    }

    virtual double dPhi(unsigned) const;

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

protected:

    struct Element_t : Detector_t::Element_t {
        Element_t(
                unsigned channel,
                unsigned adc,
                unsigned tdc
                ) :
            Detector_t::Element_t(
                channel,
                TVector3(5.1, 0, 0) // start with vector in x/y plane, is rotated in InitElements()
                ),
            ADC(adc),
            TDC(tdc)
        {}
        unsigned ADC;
        unsigned TDC;
    };

    PID(const std::vector<Element_t>& elements_init) :
        Detector_t(Detector_t::Type_t::PID),
        phi_offset0_degrees(15.476),
        elements(elements_init)
    {
        InitElements();
    }

private:
    void InitElements();
    double phi_offset0_degrees; // the offset in degrees of the first element, see InitElements()
    std::vector<Element_t> elements;

};

struct PID_2014 : PID {
    PID_2014() : PID(elements_init) {}
    virtual bool Matches(const THeaderInfo& headerInfo) const override;
    static const std::vector<Element_t> elements_init;
};

struct PID_2009_07 : PID {
    PID_2009_07() : PID(elements_init) {}
    virtual bool Matches(const THeaderInfo& headerInfo) const override;
    static const std::vector<Element_t> elements_init;
};

struct PID_2009_06 : PID {
    PID_2009_06() : PID(elements_init) {}
    virtual bool Matches(const THeaderInfo& headerInfo) const override;
    static const std::vector<Element_t> elements_init;
};

struct PID_2009_05 : PID {
    PID_2009_05() : PID(elements_init) {}
    virtual bool Matches(const THeaderInfo& headerInfo) const override;
    static const std::vector<Element_t> elements_init;
};

struct PID_2004 : PID {
    PID_2004() : PID(elements_init) {}
    virtual bool Matches(const THeaderInfo& headerInfo) const override;
    static const std::vector<Element_t> elements_init;
};


}}} // namespace ant::expconfig::detector
