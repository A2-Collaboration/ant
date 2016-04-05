#pragma once

class TH2;
class TH3;

namespace ant {
namespace MC {

struct ThetaCBToyMC {
    static double dTheta(const double theta, const double z, const double r, const double gap);
    static TH2* SimTheta(const unsigned n, const double r=25.0, const double l=10.0, const double gap=0.0, TH2* hist=nullptr);

    static void Analyse3(TH3* hist);
    static void SimThetaMulti(const unsigned n=1E+9, const double rmin=25.0, const double rmax=65, const unsigned steps=4, const double l=10.0);
};

}
}
