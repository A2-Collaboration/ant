#include "A2GeoAcceptance.h"

using namespace ant;

A2SimpleGeometry::A2SimpleGeometry():
    cb_theta_region(23.0*TMath::DegToRad(), 155.0*TMath::DegToRad()),
    cb_phi_hem1(1.5*TMath::DegToRad(),(180.0-1.5)*TMath::DegToRad()),
    taps_region(2.0*TMath::DegToRad(),19.0*TMath::DegToRad())
{

}

Detector_t::Any_t A2SimpleGeometry::DetectorFromAngles(const radian_t theta, const radian_t phi) const
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
