#pragma once
#include <vector>
#include "base/vec/vec2.h"
#include "base/math_functions/Linear.h"


namespace TPCSim {

using track_t = ant::math::LineFct;

struct diffusion_t {

    double D;

    diffusion_t(const double d):
        D(d)
    {}

    double widthAt(const double z) const {
        return D*sqrt(z);
    }

};

struct driftproperties {
    diffusion_t transverse;
    diffusion_t longitudinal;
    driftproperties(const double dt, const double dl):
        transverse(dt),
        longitudinal(dl)
    {}

};

struct tpcproperties {
    double rin=7;
    double rout=140;
    int nRings=10;
};

/**
 * @brief generate Points along a track, in the r-z plane
 * @param z0
 * @param theta
 * @param prop
 * @return
 */
std::vector<ant::vec2> generatePoints(const double z0, const double theta,
                                 const driftproperties& prop, const tpcproperties& tpc);

///@todo: add   Fit(const vector<vec2>& points)...
}



