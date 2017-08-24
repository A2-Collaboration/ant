#include "TPCSim_tools.h"
#include "TRandom2.h"

#include <iostream>

using namespace ant;
using namespace std;

namespace TPCSim {

vector<ant::vec2> generatePoints(const double z0, const double theta,
                                         const resolution_t &prop,
                                         const tpcproperties& tpc)
{
    const track_t t({0,z0},{cos(90.0-theta),sin(90.0-theta)});

    vector<vec2> points;
    points.reserve(size_t(tpc.nRings));

    const auto ringwidth = (tpc.rout-tpc.rin) / tpc.nRings;
    for(int ring = 0; ring < tpc.nRings; ++ring) {
        const double r = tpc.rin + (ring+0.5)*ringwidth;
        points.emplace_back(generatePoint(r,t,prop));
    }

    return points;
}

vec2 generatePoint(const double r, const track_t &track, const resolution_t &prop)
{
    const auto z = track(r);
    return {r, gRandom->Gaus(z, prop.dl)};
}

}

