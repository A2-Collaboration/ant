#pragma once

#include "base/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"

namespace ant {
namespace expconfig {
namespace detector {
struct PID :
        Detector_t,
        UnpackerAcquConfig // PID knows how to be filled from Acqu data
{


    virtual vec3 GetPosition(unsigned channel) const override {
        return elements[channel].Position;
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

    /**
     * @brief Get the phi angle covered by one PID element
     * @return angle in radian
     */
    virtual double dPhi(unsigned) const;

    /**
     * @brief Set an absolute phi rotation angle
     * @param offset_degrees angle in degrees
     * @todo Change to use radian
     */
    virtual void SetPhiOffset(double offset_degrees);

    /**
     * @brief Get the current phi rotation angle
     * @return angle in degrees
     * @todo Change to return values in radian
     */
    virtual double GetPhiOffset() const;

    /**
     * @brief Apply a rotation in addition to the already existing phi offset
     * @param offset The angle in degrees
     * @todo Change offset to be in radian
     */
    virtual void RotateRelative(const double offset);

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
                vec3(5.1, 0, 0) // start with vector in x/y plane, is rotated in RotateElements()
                ),
            ADC(adc),
            TDC(tdc)
        {}
        unsigned ADC;
        unsigned TDC;
    };

    /**
     * @brief Constructor
     * @param elements_init
     *
     * The default PID offset angle is conposed of 15.476 degrees from Acqu, and 7.5 degrees empirical from MC studies.
     *
     */
    PID(const std::vector<Element_t>& elements_init) :
        Detector_t(Detector_t::Type_t::PID),
        phi_offset0_degrees(15.476 + 7.5),
        elements(elements_init)
    {
        InitElements();
    }

private:
    void InitElements();
    void RotateElements();
    double phi_offset0_degrees; // the offset in degrees of the first element, see InitElements()
    std::vector<Element_t> elements;

};

struct PID_2014 : PID {
    PID_2014() : PID(elements_init) {}
    static const std::vector<Element_t> elements_init;
};

struct PID_2009_07 : PID {
    PID_2009_07() : PID(elements_init) {}
    static const std::vector<Element_t> elements_init;
};

struct PID_2009_06 : PID {
    PID_2009_06() : PID(elements_init) {}
    static const std::vector<Element_t> elements_init;
};

struct PID_2009_05 : PID {
    PID_2009_05() : PID(elements_init) {}
    static const std::vector<Element_t> elements_init;
};

struct PID_2007 : PID {
    PID_2007() : PID(elements_init) {}
    static const std::vector<Element_t> elements_init;
};

struct PID_2004 : PID {
    PID_2004() : PID(elements_init) {}
    static const std::vector<Element_t> elements_init;
};




}}} // namespace ant::expconfig::detector
