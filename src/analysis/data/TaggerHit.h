#pragma once

#include <list>
#include <memory>

#include "base/types.h"
#include "base/printable.h"
#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace data {

struct TaggerHit : printable_traits {
    index_t Channel;
    mev_t   PhotonEnergy;
    ns_t    Time;

    TaggerHit(const index_t& channel=0, const mev_t& photon_energy=0.0, const ns_t& time=0.0):
        Channel(channel),
        PhotonEnergy(photon_energy),
        Time(time) {}

    TLorentzVector GetPhotonBeam() const { return TLorentzVector(0.0, 0.0, PhotonEnergy, PhotonEnergy); }

    virtual ~TaggerHit() = default;
    std::ostream &Print(std::ostream &stream) const;
};

using TaggerHitList = std::vector<TaggerHit>;

}
}
}
