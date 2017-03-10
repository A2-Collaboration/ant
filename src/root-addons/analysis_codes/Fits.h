#pragma once

class TH1;
class TH2;
class TF1;

namespace ant {
class Fits {

public:

    struct FitResult {
        double pos;
        double N;
        double width;
        double chi2dof;
        FitResult(double pos_, double N_, double width_, double chi2dof_):
            pos(pos_),
            N(N_),
            width(width_),
            chi2dof(chi2dof_) {}
    };

    static FitResult FitEtaCalib(TH1* h, const double r_min=450.0, const double r_max=650.0);
    static FitResult FitPi0Calib(TH1* h, const double r_min=70.0, const double r_max=220.0);

    static FitResult FitPeakPol4(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max);

    static void FitSlicesPi0(TH2* h);
};

}
