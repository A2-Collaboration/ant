#pragma once
#include <string>

class TH3;
class TH2;
class TH1;
class TGraph;
class TDirectory;
class TCanvas;

namespace ant {

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

    static TH2* AnalyseChannelE(TH3* h2);

    static TH2* RelativeDiff(const TH2* h1, const TH2* h2);

    static TCanvas* getInspectorCanvas(TH2* h, const std::string& hist_base, TDirectory* dir=nullptr, const std::string& n1="");
    static TCanvas* getInspectorCanvas(TH2* h, const std::string& hist_base, TDirectory* dir1, const std::string& n1, TDirectory* dir2, const std::string& n2);
};

}
