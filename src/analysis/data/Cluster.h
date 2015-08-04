#pragma once

#include "analysis/Detector.h"

#include <limits>

namespace ant {

class Cluster {
public:
    double Energy = std::numeric_limits<double>::quiet_NaN();
    double Time = std::numeric_limits<double>::quiet_NaN();
    detector_t Detector = detector_t::None;
    unsigned CentralElement = 0;
    double RawADC = std::numeric_limits<double>::quiet_NaN();


    Cluster(double E, double t, ant::detector_t d, unsigned ch,
            double rawADC = std::numeric_limits<double>::quiet_NaN()):
        Energy(E),
        Time(t),
        Detector(d),
        CentralElement(ch),
        RawADC(rawADC)
    {}

};

}
