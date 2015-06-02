#include "A2GeoAcceptance.h"

ant::A2SimpleGeometry::A2SimpleGeometry():
    cb_theta_region(23.0*TMath::DegToRad(), 155.0*TMath::DegToRad()),
    cb_phi_hem1(1.5*TMath::DegToRad(),(180.0-1.5)*TMath::DegToRad()),
    taps_region(2.0*TMath::DegToRad(),19.0*TMath::DegToRad())
{

}

ant::detector_t ant::A2SimpleGeometry::DetectorFromAngles(const ant::radian_t theta, const ant::radian_t phi) const
{
    if( cb_theta_region.Contains(theta) ) {

        if(cb_phi_hem1.Contains(fabs(phi))) {
            return detector_t::NaI;
        }

    } else if( taps_region.Contains(theta)) {
        return detector_t::BaF2 | detector_t::PbWO4;
    }

    return detector_t::None;

}
