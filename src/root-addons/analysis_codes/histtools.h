#pragma once

#include <utility>

class TH1;
class TH1D;
class TH2;
class TGraph;

namespace ant {

struct histtools {

    /**
     * @brief Plot Integral(h1,x) / Integral(h2,x) over x
     * @param h1
     * @param h2
     * @return new histogram
     */
    static TH1* CutPosPlot(const TH1* const h1, const TH1* const h2);

    /**
     * @brief Plot Ns / sqrt (Ns / Nb)
     * @param h1
     * @param h2
     * @return new histogram
     */
    static TH1* PlotSignificance(const TH1* const h1, const TH1* const h2);

    static std::pair<TGraph*, TGraph*> MeanRMS(const TH2* const h);
    static void PlotMeansRMS(const TH2* const h);


    TH1D* ProjectX(TH2* hist);
};

}
