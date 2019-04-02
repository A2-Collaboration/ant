#include "PID.h"

#include "detail/PID_2004_elements.h"
#include "detail/PID_2007_elements.h"
#include "detail/PID_2009_05_elements.h"
#include "detail/PID_2009_06_elements.h"
#include "detail/PID_2009_07_elements.h"
#include "detail/PID_2014_elements.h"


#include "tree/TID.h"
#include "base/std_ext/math.h"

#include <cassert>
#include <cmath> // provides M_PI



using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;



double PID::dPhi(unsigned) const {
    return 2 * M_PI / elements.size();
}

void PID::SetPhiOffset(double offset_degrees) {
    phi_offset0_degrees = offset_degrees;
    RotateElements();
}

double PID::GetPhiOffset() const {
    return phi_offset0_degrees;
}

void PID::RotateRelative(const double offset) {
    SetPhiOffset(GetPhiOffset() + offset);
}

void PID::BuildMappings(vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                        vector<UnpackerAcquConfig::scaler_mapping_t>&) const
{
    for(const Element_t& element : elements)
    {
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  element.ADC);

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TDC);
    }
}

void PID::InitElements()
{
    assert(elements.size() == 24);
    for(size_t i=0; i<elements.size();i++) {
        Element_t& element = elements[i];
        // we check the order here, instead of simply computing them
        // makes handling the array mappings in detail/PID_*_elements.h easier
        if(element.Channel != i)
            throw Exception("PID element channels not in correct order");
    }
    RotateElements();
}

void PID::RotateElements()
{
    for(size_t i=0; i<elements.size();i++) {
        Element_t& element = elements[i];
        // we assume that the PID elements are ordered in positive mathematical rotation
        // the element is already initialized as a unit vector
        element.Position.SetPhi(std_ext::degree_to_radian(phi_offset0_degrees) + i*dPhi(element.Channel));
    }
}
