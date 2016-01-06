#pragma once

#include "tree/TCluster.h"

//#include "base/Detector_t.h"

//#include <list>
//#include "TVector3.h"
//#include "base/enumfield.h"

namespace ant {
namespace analysis {
namespace data {

//struct Cluster {
//    double Energy      = 0.0;
//    double ShortEnergy = 0.0;
//    double Time        = 0.0;
//    TVector3 Position;
//    Detector_t::Any_t Detector = Detector_t::Any_t::None;
//    unsigned CentralElement = 0;

//    struct Hit {
//        struct Datum {
//            Channel_t::Type_t Type;
//            double Value;
//            Datum(Channel_t::Type_t type, double value) :
//                Type(type), Value(value) {}
//        };
//        unsigned Channel;
//        std::list<Datum> Data;
//    };

//    std::list<Hit> Hits;

//    enum class Flag {
//        TouchesHole,
//        Split
//    };

//    ant::enumfield<Flag, std::int32_t> flags;

//    Cluster(double E, double ShortE, double t, Detector_t::Any_t d, unsigned ch, TVector3 position):
//        Energy(E),
//        ShortEnergy(ShortE),
//        Time(t),
//        Position(position),
//        Detector(d),
//        CentralElement(ch),
//        Hits()
//    {}


//    double GetPSARadius() const;
//    double GetPSAAngle() const;

//};

using ClusterList = std::vector<TCluster>;

}
}
}
