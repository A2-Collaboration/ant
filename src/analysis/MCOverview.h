#ifndef MCOVERVIEW_H
#define MCOVERVIEW_H

#include "AntPhysics.h"
#include "plot/Histogram.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {

class MCOverview: public ant::Physics {

    PlotList<ParticlePtr> mc_particle_stats;

    // Physics interface
public:
    MCOverview(const mev_t energy_scale=1000.0);
    virtual ~MCOverview();
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};
}
}

#endif
