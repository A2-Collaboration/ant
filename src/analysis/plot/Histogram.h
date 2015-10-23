#pragma once

#include "TColor.h"

#include <vector>
#include <string>
#include "base/interval.h"

class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {

class BinSettings: public interval<double>  {
protected:
    unsigned int bins;
public:
    BinSettings(unsigned int number_of_bins, double minimum, double maximum): interval<double>(minimum,maximum),
        bins(number_of_bins) {}

    BinSettings(unsigned int number_of_bins): interval<double>(0,number_of_bins),
        bins(number_of_bins) {}

    BinSettings(unsigned int number_of_bins, const interval<double>& i): interval<double>(i), bins(number_of_bins) {}

    virtual ~BinSettings() {}

    const unsigned int& Bins()    const { return bins; }
    unsigned int& Bins()   { return bins; }
    double BinWidth() const { return Length() / bins; }
};

class HistogramFactory {
private:
    static std::vector<EColor> colors;          // list of colors to use for 1D hitofram lines

    std::vector<EColor>::const_iterator color;  // loops over the colors to assign different ones to histograms
    unsigned int histnum;                       // number of unnamed histograms generated

    EColor GetNextColor();
    bool loopColors;
    unsigned int GetNextHistnum();


    static HistogramFactory* instance;

public:
    std::string GetNextHistName(const std::string &name="");

    HistogramFactory();

    /**
     * @brief Make a new 1D histogram
     * @param title The title of the histogram
     * @param xlabel X axis label (required, cause... always label your axes ;-) )
     * @param ylabel Y axis label (required)
     * @param bins Binning settings, (Number, min, max)
     * @param name optional name. If "", a number will be inserted
     * @return A pointer to the new histogram
     */
    template<class Hist = TH1D>
    Hist* Make1D (const std::string& title, const std::string& xlabel, const std::string& ylabel, const BinSettings& bins=BinSettings(100,0,100), const std::string& name="")
    {

        Hist* h = new Hist( GetNextHistName(name).c_str(), title.c_str(), bins.Bins(), bins.Start(), bins.Stop());
        h->SetXTitle(xlabel.c_str());
        h->SetYTitle(ylabel.c_str());
        h->SetLineColor(GetNextColor());
        return h;
    }

    /**
     * @brief Make a new 2D histogram
     * @param title The title of the histogram
     * @param xlabel X axis label (required, cause... always label your axes ;-) )
     * @param ylabel Y axis label (required)
     * @param xbins x axis binning settings, (Number, min, max)
     * @param ybins y axis binning settings, (Number, min, max)
     * @param name optional name. If "", a number will be inserted
     * @return A pointer to the new histogram
     */
    TH2D* Make2D (const std::string& title,
                  const std::string& xlabel,
                  const std::string& ylabel,
                  const BinSettings& xbins=BinSettings(100),
                  const BinSettings& ybins=BinSettings(100),
                  const std::string& name="");

    TH3D* Make3D (const std::string& title,
                  const std::string& xlabel,
                  const std::string& ylabel,
                  const std::string& zlabel,
                  const BinSettings& xbins=BinSettings(100),
                  const BinSettings& ybins=BinSettings(100),
                  const BinSettings& zbins=BinSettings(100),
                  const std::string& name="");

    void ApplySettings(TH1D* hist, const std::string& title="", const std::string& xlabel="", const std::string& ylabel="");
    void ApplySettings(TH2D* hist, const std::string& title="", const std::string& xlabel="", const std::string& ylabel="");

    void SetLoopColors( bool onoff ) { loopColors=onoff; }
    bool GetLoopColors() { return loopColors; }

    void ResetColors();
    void SetOutputRoot(TDirectory* dir);

    static HistogramFactory& Default() { return *instance; }

};

}
}
