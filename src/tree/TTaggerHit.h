#pragma once

#include "TID.h" // TKeyValue

#include "base/vec/LorentzVec.h"

#include <vector>
#include <iomanip>
#include <sstream>

namespace ant {

struct TTaggerHit
{
    unsigned Channel;
    double   PhotonEnergy;
    double   Time;

    struct Electron_t {
        unsigned Channel;
        double Timing;
        double QDCEnergy;

        Electron_t(unsigned channel, double timing,
                   double energy = std::numeric_limits<double>::quiet_NaN()) :
            Channel(channel), Timing(timing), QDCEnergy(energy)
        {}

        Electron_t() = default;
        template<class Archive>
        void serialize(Archive& archive) {
            archive(Channel, Timing, QDCEnergy);
        }

        friend std::ostream& operator<<( std::ostream& s, const Electron_t& o) {
            s << "(" << o.Timing;
            if(std::isfinite(o.QDCEnergy))
                s << "," << o.QDCEnergy;
            return s << ")";
        }
    };

    std::vector<Electron_t> Electrons;

    TTaggerHit(unsigned channel, double photonE, double time,
               double qdc_energy = std::numeric_limits<double>::quiet_NaN()) :
        Channel(channel),
        PhotonEnergy(photonE),
        Time(time),
        Electrons{ {channel, time, qdc_energy}  }
    {}

    TTaggerHit() = default;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(PhotonEnergy, Time, Channel, Electrons);
    }

    LorentzVec GetPhotonBeam() const {
        return {{0.0, 0.0, PhotonEnergy}, PhotonEnergy};
    }


    friend std::ostream& operator<<( std::ostream& s, const TTaggerHit& o) {
        return s << "TTaggerHit: Electrons=" << o.Electrons.size()
                 << " Channel=" << o.Channel << " PhotonEnergy=" << o.PhotonEnergy
                 << " Time=" << o.Time;
    }


};

}
