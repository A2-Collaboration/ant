#pragma once

#include "TID.h" // TKeyValue

#include "base/printable.h"

#include <vector>
#include <iomanip>
#include <sstream>

namespace ant {

struct TTaggerHit : printable_traits

{
    double PhotonEnergy;
    double Time;
    unsigned Channel;
    std::vector< TKeyValue<double> > Electrons; // Key=channel, Value=timing


    TTaggerHit(double photonE, const TKeyValue<double>& hit) :
        PhotonEnergy(photonE),
        Time(hit.Value),
        Channel(hit.Key),
        Electrons{hit}
    {}
    TTaggerHit(unsigned channel, double photonE, double time) :
        PhotonEnergy(photonE),
        Time(time),
        Channel(channel)
    {}

    TTaggerHit() {}
    virtual ~TTaggerHit() {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(PhotonEnergy, Time, Channel, Electrons);
    }


    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TTaggerHit: Electrons=" << Electrons.size()
                 << " PhotonEnergy=" << PhotonEnergy << " Time=" << Time;
    }


};

}
