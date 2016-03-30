#pragma once

class TH1;


namespace ant {

struct histtools {

    /**
     * @brief Plot Integral(h1,x) / Integral(h2,x) over x
     * @param h1
     * @param h2
     * @return new histogram
     */
    static TH1* CutPosPlot(const TH1* const h1, const TH1* const h2);
};

}
