#pragma once
#include <Rtypes.h>

class TH1;
class TH1D;
class TH2;
class TH2D;
class TF1;
class TCanvas;

namespace ant {
struct hstack;

class Omega {

public:

    struct FitResult {
        double pos;
        double N;
        double width;
        double chi2dof;
        FitResult(double pos_=0., double N_=0., double width_=0., double chi2dof_=-1.):
            pos(pos_),
            N(N_),
            width(width_),
            chi2dof(chi2dof_) {}
    };

    enum FCT_t {
        eGAUS,
        eCrystalBall
    };

    static FitResult FitHist(TH1* h, const double omega_mass_ = -1.0, const bool fixOmegaMass=false, const FCT_t fct = eGAUS, const double r_min=650.0, const double r_max=900.0);
    static FitResult FitHistCrystalBall(TH1* h, const bool fixOmegaMass=false, const double r_min=650.0, const double r_max=900.0);

//    static TF1* nCB();
    static TF1* nCB_pol3();

    static TF1* AsymGaus_pol3();

    static void FitOmegaPeak(const bool fixOmegaMass=false, const double r_min=650.0, const double r_max=900.0);

    static void FitOmega2D(TH2* h);

    static TH2D* SampleDiffXsectionOmega();
    static TH2D* SampleDiffXsectionEta();
    static TH2D* SampleDiffXsectionPi0();

    static TH1D* getSignalYield(const std::string& meson = "eta");

    static void SaveStacks(const std::string& path_spec, const std::string &fname, const int start=0, const int stop=9);

    static void PlotFitted(const std::string& file);

    static TH2D* ConvertMesonCountHistToCrossSec(const TH2D* mesonCounts);

    static TCanvas* Printhstack(const std::string& name);
    static TCanvas* Printhstack(ant::hstack* hs);

    static void RoundBins(TH2* h, int d);

    static TH1D* HistOfBins(const TH2* h2);

    static std::string TH1ToLaTeX(const TH1* h, const int p=3);
};

}
