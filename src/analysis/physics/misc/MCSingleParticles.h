#ifndef MCSINGLEPARTICLES_H
#define MCSINGLEPARTICLES_H

#include "AntPhysics.h"
#include "plot/Histogram.h"


namespace ant {
namespace analysis {

class MCSingleParticles: public ant::Physics {

public:
    typedef std::pair<const TrackPtr&, const ParticlePtr&> Track_MC_pair;
    typedef std::pair<const TrackList&, const ParticlePtr&> MC_tracklist_pair;
    typedef std::pair<const ParticlePtr&, const ParticlePtr&> Rec_MC_pair;
private:

    PlotList<Track_MC_pair> MC_track_pair_stats;
    PlotList<MC_tracklist_pair> MC_tracklist_pair_stats;
    PlotList<Rec_MC_pair> Rec_MC_stats;

    // Physics interface
public:
    MCSingleParticles(const mev_t energy_scale=1000.0);
    virtual ~MCSingleParticles();
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};
}
}

#endif
