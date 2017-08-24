#pragma once
#include <vector>
#include "base/vec/vec2.h"
#include <ostream>
#include "base/math_functions/Linear.h"

#include "APLCON.hpp"

namespace ant {

namespace TPCSim {

using track_t = ant::math::LineFct;
struct Value_t {

    constexpr Value_t(double v, double s) : V_S_P{v, s, 0}{}

    std::tuple<double,double,double> V_S_P; // short for Value, Sigma, Pull

    bool Fixed = false;
    bool Poisson = false;

    // fitter is able to access fields to fit them (last one pull is just written, actually)
    template<size_t N>
    std::tuple<double&> linkFitter() noexcept {
        return std::tie(std::get<N>(V_S_P));
    }

    // user can provide settings for each variabl
    template<size_t innerIdx>
    APLCON::Variable_Settings_t getFitterSettings(size_t outerIdx) const noexcept {
        (void)outerIdx; // unused, provided to user method for completeness
        APLCON::Variable_Settings_t settings;
        if(Fixed)
            settings.StepSize = 0;
        if(Poisson)
            settings.Distribution = APLCON::Distribution_t::Poissonian;
        return settings;
    }

    // the following methods are just user convenience (not required by fitter)

    friend std::ostream& operator<<(std::ostream& s, const Value_t& o) {
        return s << "(" << o.Value() << "," << o.Sigma() << "," << o.Pull() << ")";
    }

    double Value() const {
        return std::get<APLCON::ValueIdx>(V_S_P);
    }

    double Sigma() const {
        return std::get<APLCON::SigmaIdx>(V_S_P);
    }

    double Pull() const {
        return std::get<APLCON::PullIdx>(V_S_P);
    }

    operator double() const {
        return std::get<APLCON::ValueIdx>(V_S_P);
    }
};


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
