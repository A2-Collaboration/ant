#pragma once

class TH3;
class TH2;
class TH1;

namespace ant {

struct PeakFitResult_t {
    double chi2dof;
    double pos;
    double sigma;

    PeakFitResult_t(const double chi2dof_=0.0, const double pos_=0.0, const double sigma_=0.0):
        chi2dof(chi2dof_), pos(pos_), sigma(sigma_) {}
};

struct TowPi0_MCSmearing_Tool {

    static PeakFitResult_t Fit(TH1* hist);

};

}
