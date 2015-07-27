#pragma once

#include "analysis/Detector.h"

namespace ant {

class Cluster {
public:
    double Energy = 0.0;
    double Time = 0.0;
    detector_t Detector = detector_t::None;
    unsigned CentralElement = 0;

    Cluster(double E, double t, ant::detector_t d, unsigned ch):
        Energy(E),
        Time(t),
        Detector(d),
        CentralElement(ch) {}

};

}
