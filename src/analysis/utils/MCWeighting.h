#pragma once

#include "plot/HistogramFactories.h"
#include "tree/TParticle.h"
#include "base/WrapTTree.h"

#include <vector>
#include <map>

namespace ant {
namespace analysis {
namespace utils {

class MCWeighting {

protected:
    struct coefficients_t {
        interval<double>    BeamE;
        std::vector<double> LegendreCoefficients;

        coefficients_t(const interval<double>& beamE, const std::vector<double> coeffs) :
            BeamE(beamE), LegendreCoefficients(coeffs)
        {}

        bool operator<(const coefficients_t& o) const {
            return BeamE.Center() < o.BeamE.Center();
        }
    };
    using database_t = std::vector<coefficients_t>;

    static database_t SanitizeDatabase(database_t d);

    double GetBeamE(const TParticleTree_t& tree);
    double GetCosTheta(const TParticleTree_t& tree);
    double GetN(double beamE, double cosTheta) const;

    const database_t& Database;

    unsigned nParticleTrees = 0;
    double N_sum = 0;
    double last_N = std_ext::NaN;


public:

    static const database_t EtaPrime;

    MCWeighting(const HistogramFactory& HistFac, const database_t& database);

    void SetParticleTree(const TParticleTree_t& tree);



    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };
};

}
}
}