#pragma once

#include "base/Detector_t.h"

#include <list>
#include "TVector3.h"
#include "base/enumfield.h"

namespace ant {
namespace analysis {
namespace data {

class Cluster {
public:
    double Energy = 0.0;
    double Time = 0.0;
    TVector3 pos;
    Detector_t::Any_t Detector = Detector_t::Any_t::None;
    unsigned CentralElement = 0;

    struct Hit {
        struct Datum {
            Channel_t::Type_t Type;
            double Value;
            Datum(Channel_t::Type_t type, double value) :
                Type(type), Value(value) {}
        };
        unsigned Channel;
        std::list<Datum> Data;
    };

    std::list<Hit> Hits;

    enum class Flag {
        TouchesHole,
        Split
    };

    ant::enumfield<Flag, std::int32_t> flags;

    Cluster(double E, double t, Detector_t::Any_t d, unsigned ch, TVector3 position):
        Energy(E),
        Time(t),
        pos(position),
        Detector(d),
        CentralElement(ch),
        Hits()
    {}

};

using ClusterList = std::vector<Cluster>;

}
}
}
