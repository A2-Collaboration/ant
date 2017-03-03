#pragma once

#include "base/interval.h"
#include "base/Detector_t.h"

namespace ant {
namespace analysis {
namespace utils {


class A2SimpleGeometry {
protected:
    const IntervalD cb_theta_region;  // theta
    const IntervalD cb_phi_hem1;      // phi region for one hemisphere
    const IntervalD taps_region;      // theta

public:
    A2SimpleGeometry();

    Detector_t::Any_t DetectorFromAngles( const radian_t theta, const radian_t phi) const;

    template <class T>
    Detector_t::Any_t  DetectorFromAngles( const T& v) const {
        return DetectorFromAngles(v.Theta(), v.Phi());
    }

};

}}} // namespace ant::analysis::utils
