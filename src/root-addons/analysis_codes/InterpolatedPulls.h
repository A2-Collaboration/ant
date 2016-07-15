#pragma once
#include <string>
#include <list>
#include "TDirectory.h"

class TH2;
class DirList;
class TCanvas;

namespace ant {

struct InterpolatedPulls {

    static void PlotDirectory(std::list<TDirectory*> dirs, const std::string& prefix, const std::string& title="");
    static void PlotComparePulls2();
    static void PlotPullSigmas(const std::string& treename, TDirectory* dir);

    static void PlotComparePulls(TDirectory* red, TDirectory* blue);
    static void PlotAllSigams(TDirectory* dir=gDirectory);

    static void TestInterpolation(const std::string& filename);

    static TCanvas* ConvergencePlots(std::list<TH2*> hists);
    static void ConvergencePlots2();

};

struct ConvergencePlot {

    DirList* list=nullptr;

    ConvergencePlot();
    ~ConvergencePlot();

    void Add(TDirectory* dir);

    void Plot(const double min=0.0, const double max=0.0, const int ww=-1, const int wh=-1);

};

}
