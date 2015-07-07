#ifndef TAGGERHIT_H
#define TAGGERHIT_H

#include <list>
#include <memory>

#include "base/types.h"
#include "base/printable.h"
#include "TLorentzVector.h"

namespace ant {

class TaggerHit: public ant::printable_traits {
private:
    index_t channel;
    mev_t   photon_energy;
    ns_t  time;

public:
    TaggerHit(const index_t& _channel=0, const mev_t& _photon_energy=0.0, const ns_t& _time=0.0):
        channel(_channel),
        photon_energy(_photon_energy),
        time(_time) {}

    virtual ~TaggerHit() {}

    index_t& Channel() { return channel; }
    const index_t& Channel() const { return channel; }

    mev_t& PhotonEnergy() { return photon_energy; }
    const mev_t& PhotonEnergy() const { return photon_energy; }

    ns_t& Time() { return time; }
    const ns_t& Time() const { return time; }

    TLorentzVector PhotonBeam() const { return std::move(TLorentzVector(0.0, 0.0, PhotonEnergy(), PhotonEnergy())); }


    std::ostream &Print(std::ostream &stream) const;
};

using TaggerHitPtr   = std::shared_ptr<TaggerHit>;
using TaggerHistList = std::vector<TaggerHitPtr>;

}

#endif
