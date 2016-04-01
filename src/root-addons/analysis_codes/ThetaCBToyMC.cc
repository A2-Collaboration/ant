#include "ThetaCBToyMC.h"
#include "base/std_ext/math.h"
#include "TH2.h"
#include "TH2D.h"


using namespace ant;
using namespace ant::MC;
using namespace ant::std_ext;

#include <cmath>
#include "TRandom2.h"

double ThetaCBToyMC::dTheta(const double theta, const double z, const double r) {
    const auto a0 = r*cos(theta);
    const auto b  = r*sin(theta);

    return atan2(b,a0-z);
}

TH2* ThetaCBToyMC::SimTheta(const unsigned n, const double r, const double l, TH2* hist)
{
    const double theta_min = degree_to_radian(20.0);
    const double theta_max = degree_to_radian(160.0);

    if(!hist) {
        hist = new TH2D("thetacb", Form("d#theta: %f cm target, %f cm radius", l, r), 180, radian_to_degree(theta_min), radian_to_degree(theta_max), 200,-20, 20);
        hist->SetXTitle("#theta [#circ]");
        hist->SetYTitle("d#theta [#circ]");
    }

    TRandom2 rng;
    rng.SetSeed();

    for(unsigned i=0; i<n; ++i) {
        const auto theta = rng.Uniform(theta_max - theta_min) + theta_min;
        const auto z     = rng.Uniform(l) - l/2.0;

        const auto dtheta = theta - dTheta(theta, z, r);

        hist->Fill(radian_to_degree(theta), radian_to_degree(dtheta));
    }

    return hist;

}
