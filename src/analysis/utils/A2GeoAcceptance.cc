#include "A2GeoAcceptance.h"

#include "base/std_ext/math.h"

using namespace ant;
using namespace ant::analysis::utils;

A2SimpleGeometry::A2SimpleGeometry():
    cb_theta_region(std_ext::degree_to_radian(IntervalD{23.0, 155.0})),
    cb_phi_hem1(std_ext::degree_to_radian    (IntervalD{1.5,(180.0-1.5)})),
    taps_region(std_ext::degree_to_radian    (IntervalD{2.0,19.0}))
{

}

Detector_t::Any_t A2SimpleGeometry::DetectorFromAngles(const radian_t theta,
                                                       const radian_t phi) const
{
    if( cb_theta_region.Contains(theta) ) {

        if(cb_phi_hem1.Contains(fabs(phi))) {
            return Detector_t::Type_t::CB;
        }

    } else if( taps_region.Contains(theta)) {
        return Detector_t::Type_t::TAPS;
    }

    return Detector_t::Any_t::None;

}
