#pragma once
#include <string>

class TH3;
class TH2;
class TH1;
class TGraph;
class TDirectory;

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
    TGraph* g;
    TH1* h;
    MultiChannelFitResult_t(TGraph* graph, TH1* hist):g(graph), h(hist) {}
};

struct TowPi0_MCSmearing_Tool {

    static PeakFitResult_t Fit(TH1* hist, const std::string& prefix, const bool verbose=false);

    static void DrawSame(TH1* h1, TH1* h2);

    static channelFitResult_t FitChannel(TH2* h2, const int ch);

    static MultiChannelFitResult_t FitAllChannels(TH2* h2, const std::string& prefix);

    static void CompareMCData(TDirectory* mc, TDirectory* data);

    static TH2* AnalyseChannelE(TH3* h2);

    static TH2* RelativeDiff(const TH2* h1, const TH2* h2);
};

}
