#pragma once

#include "base/Detector_t.h"

#include <list>

namespace ant {


class Cluster {
public:
    double Energy = 0.0;
    double Time = 0.0;
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

    Cluster(double E, double t, Detector_t::Any_t d, unsigned ch):
        Energy(E),
        Time(t),
        Detector(d),
        CentralElement(ch),
        Hits()
    {}

};

using ClusterList = std::vector<Cluster>;

}
