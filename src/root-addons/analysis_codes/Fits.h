#pragma once

class TH1;
class TH2;
class TH3;
class TF1;

namespace ant {
class Fits {

public:

    struct FitResult {
        double pos;
        double N;
        double width;
        double chi2dof;
        TF1* sum = 0;
        TF1* bkg = 0;
        TF1* sig = 0;
        double position_error;
        FitResult(double pos_, double N_, double width_, double chi2dof_, TF1* Sum=0, TF1* Bkg=0, TF1* Sig=0, double pos_error = 0):
            pos(pos_),
            N(N_),
            width(width_),
            chi2dof(chi2dof_), sum(Sum), bkg(Bkg),  sig(Sig), position_error(pos_error)  {}
    };

    static FitResult FitEtaCalib(TH1* h, const double r_min=450.0, const double r_max=650.0);
    static FitResult FitEtaPrimeCalib(TH1* h, const double r_min=800.0, const double r_max=1050.0);
    static FitResult FitPi0Calib(TH1* h, const double r_min=70.0, const double r_max=220.0);
    static FitResult FitPi0Calib0(TH1* h, const double r_min=70.0, const double r_max=220.0);
    static FitResult FitPi0CalibGaussian(TH1* h, const double r_min=70.0, const double r_max=220.0);

    static FitResult FitPeakPol0(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max);
    static FitResult FitPeakPol4(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max);
    static FitResult FitPeakCrystalBallPol0(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max);
    static FitResult FitPeakCrystalBallPol4(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max);
    static FitResult FitPeakCrystalBallPol6(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max);
    static FitResult FitPeakPol6(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max);

    static void FitSlicesPi0(TH2* h);
    static void FitSlicesZVertex(TH3* h);
    static void FitSlicesEta(TH2* h);
    static void FitSlicesEtaPrime(TH2* h);
    static void FitSlicesAlpha(TH2* h);

};

}
