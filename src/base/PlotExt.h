#pragma once


#include "base/interval.h"
#include <vector>

class TH1;
class TH2;
class TH2D;
class TGraph;
class TGraphErrors;


namespace ant{

double GetMaxPos(TH1* hist);

struct TH2Ext {
    static void MakeSameZRange(std::vector<TH2*> hists);
    static void MakeSameZRange(std::vector<TH2*> hists, const ant::interval<double>& range);

    static void ClearHistogram(TH2D* hist, const double v=0.0);

    static interval<double> getZRange(const TH2& hist);
};

struct GraphExt
{
    static std::size_t FillGraph(TGraph* graph, const double x, const double y);
    static std::size_t FillGraphErrors(TGraphErrors* graph,
                                       const double x, const double y,
                                       const double xerr, const double yerr);
};

}
