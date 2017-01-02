#pragma once
#include <string>
#include <list>
#include "TCanvas.h"

class TH3;
class TH2;
class TH1;
class TGraph;
class TDirectory;
class TFile;

namespace ant {

class TCalibrationData;

struct PeakFitResult_t {
    double chi2dof;
    double pos;
    double sigma;
    int status;

    PeakFitResult_t(const double chi2dof_=0.0, const double pos_=0.0, const double sigma_=0.0, const int status_=0):
        chi2dof(chi2dof_), pos(pos_), sigma(sigma_), status(status_) {}
};

struct channelFitResult_t {
    PeakFitResult_t result;
    TH1* h;
    channelFitResult_t(const PeakFitResult_t& r, TH1* hist): result(r), h(hist) {}
};

struct MultiChannelFitResult_t {
    TH1* position;
    TH1* sigma;
    MultiChannelFitResult_t(TH1* pos, TH1* sig):position(pos), sigma(sig) {}
};

struct TwoPi0_MCSmearing_Tool {

    static PeakFitResult_t Fit(TH1* hist, const std::string& prefix, const bool verbose=false);

    static void DrawSame(TH1* h1, TH1* h2);

    static channelFitResult_t FitChannel(TH2* h2, const int ch);

    static MultiChannelFitResult_t FitAllChannels(TH2* h2);

    static void CompareMCData1D(TDirectory* mc, TDirectory* data);
    static void CompareMCData2D(TDirectory* mc, TDirectory* data, const std::string& folder="ETheta");

    /**
     * @brief Calculate new energy smearing factors
     * @param sigma_data 2D histogram containing the peak widths obained from exp. data (stays fixed over iterations)
     * @param current_sigma_MC 2D histogram containing the peak widths obtained from the most recent iteration step on MC data
     * @param last_diff difference data - smeard_mc from last iteration
     * @param sum_diffs sum of the differences (data-smeard_mc) from all prev. iterations
     * @return A new histogram containing the energy smearing for the next iteration. BinContent = -1.0 means no data
     */
    static TH2* CalculateUpdatedSmearing(const TH2* sigma_data, const TH2* current_sigma_MC, const TH2* last_smear);

    /**
     * @brief Calculate Initial eneryg Smearing
     * @param sigma_data
     * @param current_sigma_MC
     * @return A new histogram containing the energy smearing for the first iteration. BinContent = -1.0 means no data
     */
    static TH2* CalculateInitialSmearing(const TH2* sigma_data, const TH2* sigma_MC);

    static TH2* AnalyseChannelE(TH3* h2);

    static TH2* RelativeDiff(const TH2* h1, const TH2* h2);

    static TCanvas* getInspectorCanvas(TH2* h, const std::string& hist_base, TDirectory* dir=nullptr, const std::string& n1="");
    static TCanvas* getInspectorCanvas(TH2* h, const std::string& hist_base, TDirectory* dir1, const std::string& n1, TDirectory* dir2, const std::string& n2);

    static TH2* Decode(const TCalibrationData& cdata);
    static TH2* LoadAndDecode(TFile* f);
};

class TBinGraphCanvas : public TCanvas {
public:

    TBinGraphCanvas();
    virtual ~TBinGraphCanvas();

    void Add(TH2* h);

    void HandleInput(const EEventType button, const Int_t px, const Int_t py) override;

private:
    std::list<TH2*> hists;
    TGraph* graph;
    void FillGraph(const int x, const int y);
    int obx;
    int oby;
};

}
