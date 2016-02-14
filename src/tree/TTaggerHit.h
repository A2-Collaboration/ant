#pragma once

#include "TID.h" // TKeyValue

#include "base/printable.h"

#include "TLorentzVector.h"

#include <vector>
#include <iomanip>
#include <sstream>

namespace ant {

struct TTaggerHit : printable_traits

{
    unsigned Channel;
    double   PhotonEnergy;
    double   Time;

    struct Electron_t : printable_traits {
        double Timing;
        double QDCEnergy;

        Electron_t(double timing,
                   double energy = std::numeric_limits<double>::quiet_NaN()) :
            Timing(timing), QDCEnergy(energy)
        {}

        Electron_t() {}
        virtual ~Electron_t() {}
        template<class Archive>
        void serialize(Archive& archive) {
            archive(Timing, QDCEnergy);
        }

        virtual std::ostream& Print( std::ostream& s) const override {
            s << "(" << Timing;
            if(std::isfinite(QDCEnergy))
                s << "," << QDCEnergy;
            return s << ")";
        }
    };

    std::vector< TKeyValue<Electron_t> > Electrons; // Key=channel

    TTaggerHit(unsigned channel, double photonE, double time,
               double qdc_energy = std::numeric_limits<double>::quiet_NaN()) :
        Channel(channel),
        PhotonEnergy(photonE),
        Time(time),
        Electrons{ {channel, {time, qdc_energy} } }
    {}

    TTaggerHit() {}
    virtual ~TTaggerHit() {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(PhotonEnergy, Time, Channel, Electrons);
    }

    TLorentzVector GetPhotonBeam() const {
        return TLorentzVector(0.0, 0.0, PhotonEnergy, PhotonEnergy);
    }


    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TTaggerHit: Electrons=" << Electrons.size()
                 << " Channel=" << Channel << " PhotonEnergy=" << PhotonEnergy
                 << " Time=" << Time;
    }


};

}
