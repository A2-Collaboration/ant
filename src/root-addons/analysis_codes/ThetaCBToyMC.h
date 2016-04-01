#pragma once

class TH2;

namespace ant {
namespace MC {

struct ThetaCBToyMC {
    static double dTheta(const double theta, const double z, const double r);
    static TH2* SimTheta(const unsigned n, const double r=25.0, const double l=10.0, TH2* hist=nullptr);
};

}
}
