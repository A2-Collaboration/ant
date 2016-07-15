#pragma once


#include "base/interval.h"
#include <vector>

class TH1;
class TH2;
class TH2D;

namespace ant{

double GetMaxPos(TH1* hist);

struct TH2Ext {
    static void MakeSameZRange(std::vector<TH2*> hists);
    static void MakeSameZRange(std::vector<TH2*> hists, const ant::interval<double>& range);

    static void ClearHistogram(TH2D* hist, const double v=0.0);

    static interval<double> getZRange(const TH2& hist);
};

}
