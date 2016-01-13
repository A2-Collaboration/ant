#ifndef DETECTOR_PLOTS_H
#define DETECTOR_PLOTS_H

#include <string>

namespace ant {
    class TH2CB;
    class TH2TAPS;
}


class DetectorPlots {
public:
    static ant::TH2TAPS* MakeTAPSTheta(const std::string& setup_name);
    static ant::TH2TAPS* MakeTAPSPhi(const std::string& setup_name);

    static ant::TH2CB*   MakeCBTheta(const std::string& setup_name);
    static ant::TH2CB*   MakeCBPhi(const std::string& setup_name);

    static void          PlotCBTheta(const std::string& setup_name);
    static void          PlotCBPhi(const std::string& setup_name);

    static void          PlotTAPSTheta(const std::string& setup_name);
    static void          PlotTAPSPhi(const std::string& setup_name);

};


#endif
