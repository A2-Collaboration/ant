#pragma once
#include <string>
#include "TDirectory.h"

class TH2;

namespace ant {

struct InterpolatedPulls {

    static void PlotDirectory(TDirectory* dir, const std::string& prefix, const int cols, const int rows, const std::string& title="", TDirectory* dir2=nullptr);
    static void PlotPullSigmas(const std::string& treename, TDirectory* dir);

    static void PlotComparePulls(TDirectory* red, TDirectory* blue);
    static void PlotAllSigams(TDirectory* dir=gDirectory);

};

}
