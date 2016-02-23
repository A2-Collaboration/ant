#pragma once

#include "physics/Physics.h"

#include "utils/A2GeoAcceptance.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class MCTrueAcceptance: public Physics {
protected:
    TH1D* detect;
    utils::A2SimpleGeometry geo;
    unsigned int events_seen;

    struct det_hit_count_t {
        det_hit_count_t(unsigned int _cb=0, unsigned int _taps=0): cb(_cb), taps(_taps) {}
        unsigned int cb;
        unsigned int taps;
    };

    det_hit_count_t AllAccepted(const TParticleList& particles);

    bool alldetectable(const TParticleList& particles) const;

public:
    MCTrueAcceptance(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
