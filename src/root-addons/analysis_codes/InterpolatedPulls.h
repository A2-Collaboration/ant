#pragma once
#include <string>
#include "TDirectory.h"

class TH2;
class DirList;

namespace ant {

struct InterpolatedPulls {

    static void PlotDirectory(TDirectory* dir, const std::string& prefix, const std::string& title="", TDirectory* dir2=nullptr);
    static void PlotPullSigmas(const std::string& treename, TDirectory* dir);

    static void PlotComparePulls(TDirectory* red, TDirectory* blue);
    static void PlotAllSigams(TDirectory* dir=gDirectory);

    static void TestInterpolation(const std::string& filename);

};

struct ConvergencePlot {

    DirList* list=nullptr;

    ConvergencePlot();
    ~ConvergencePlot();

    void Add(TDirectory* dir);

    void Plot();

};

}
