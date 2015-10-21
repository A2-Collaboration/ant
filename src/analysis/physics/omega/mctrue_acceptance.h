#pragma once

#include "physics/Physics.h"
#include "plot/SmartHist.h"

#include "utils/A2GeoAcceptance.h"

namespace ant {
namespace analysis {
namespace physics {

class MCTrueAcceptance: public Physics {
protected:
    SmartHist1<std::string> detect;
    utils::A2SimpleGeometry geo;
    unsigned int events_seen;

    struct det_hit_count_t {
        det_hit_count_t(unsigned int _cb=0, unsigned int _taps=0): cb(_cb), taps(_taps) {}
        unsigned int cb;
        unsigned int taps;
    };

    det_hit_count_t AllAccepted(const data::ParticleList& particles);

    bool alldetectable(const data::ParticleList& particles) const;

public:
    MCTrueAcceptance(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void Finish() override;
    void ShowResult() override;
};

}
}
}
