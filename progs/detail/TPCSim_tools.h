#pragma once
#include <vector>
#include "base/vec/vec2.h"

namespace ant {

namespace TPCSim {

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

std::vector<vec2> generatePoints(const double z0, const double theta,
                                 const driftproperties& prop);

}

///@todo: add   Fit(const vector<vec2>& points)...
}
