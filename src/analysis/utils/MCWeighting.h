#pragma once

#include "analysis/plot/HistogramFactory.h"
#include "tree/TParticle.h"
#include "base/WrapTTree.h"

#include <vector>
#include <map>

namespace ant {
namespace analysis {
namespace utils {

class MCWeighting {
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

public:

    struct item_t {
        const ParticleTypeDatabase::Type& Type;
        database_t Database;
    };

    static const std::string treeName;
    static const item_t EtaPrime;
    static const item_t Omega;
    static const item_t Pi0;

    MCWeighting(const HistogramFactory& histFac, const item_t& item);

    // usage of those methods is tricky...see also test/TestMCWeighting physics class
    // 1) SetParticleTree should be called for each event encountered
    // 2) Fill() should be called everytime a tree is filled with "physics" results
    // 3) Finish() must be called after no more Fill/SetParticleTree are done
    void SetParticleTree(const TParticleTree_t& tree);
    void Fill();
    void Finish();

    // can be used to access the "internal" MCWeights tree
    // alternatively, get the tree by its name MCWeighting::treeName
    bool FriendTTree(TTree* tree);

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    struct tree_t : WrapTTree {
        // use "unique" branch name to make it easy to friend
        // this TTree
        ADD_BRANCH_T(double, MCWeight)
    };

//protected:

    static database_t SanitizeDatabase(database_t d);

    double GetBeamE(const TParticleTree_t& tree);
    double GetCosTheta(const TParticleTree_t& tree);
    double GetN(double beamE, double cosTheta) const;

    const item_t& Item;

    unsigned nParticleTrees = 0;
    double N_sum = 0;
    double last_N = std_ext::NaN;

    HistogramFactory HistFac;
    tree_t t;
    TTree* treeWeighted = nullptr;
};

}
}
}
