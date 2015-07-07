#ifndef A2GEOACCEPTANCE_H
#define A2GEOACCEPTANCE_H

#include "base/interval.h"
#include "TMath.h"
#include "Detector.h"

namespace ant {

class A2SimpleGeometry {
protected:
    IntervalD cb_theta_region;  // theta
    IntervalD cb_phi_hem1;      // phi region for one hemisphere
    IntervalD taps_region;      // theta

public:
    A2SimpleGeometry();

    detector_t DetectorFromAngles( const radian_t theta, const radian_t phi) const;

    template <class T>
    detector_t DetectorFromAngles( const T& v) const {
        return DetectorFromAngles(v.Theta(), v.Phi());
    }

};

}
#endif
