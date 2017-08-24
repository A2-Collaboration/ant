#include "TPCSim_tools.h"
#include <iostream>

using namespace ant;
using namespace std;

namespace TPCSim {

vector<ant::vec2> generatePoints(const double z0, const double theta,
                                         const driftproperties&,
                                         const tpcproperties& tpc)
{
    const track_t t({0,z0},{cos(90.0-theta),sin(90.0-theta)});

    const auto ringwidth = (tpc.rout-tpc.rin) / tpc.nRings;
    for(int ring = 0; ring < tpc.nRings; ++ring) {
        const double r = tpc.rin + (ring+0.5)*ringwidth;
        cout << r << endl;
    }

    return {};
}

}

