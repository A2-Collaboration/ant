#pragma once

#include "tree/TParticle.h"

#include <vector>
#include <map>

namespace ant {
namespace analysis {
namespace utils {

class MCWeighting {

    struct coefficients_t {
        double BeamE;
        std::vector<double> LegendreCoefficients;

        coefficients_t(double beamE, const std::vector<double> coeffs) :
            BeamE(beamE), LegendreCoefficients(coeffs)
        {}

        bool operator<(const coefficients_t& o) const {
            return BeamE < o.BeamE;
        }
    };
    using database_t = std::map<const ParticleTypeDatabase::Type*, std::vector<coefficients_t>>;
    const database_t database;
    static database_t MakeDatabase();

public:

    MCWeighting();

    double GetWeight(const TParticleTree_t& tree) const;

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };
};

}
}
}