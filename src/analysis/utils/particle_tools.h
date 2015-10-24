#pragma once

#include "analysis/data/Particle.h"
#include "base/ParticleTypeTree.h"

#include <string>

#include "TLorentzVector.h"

class TH1;
class TTree;

namespace ant {
namespace analysis {
namespace utils {

struct ParticleVars {
    double E;
    double Theta;
    double Phi;
    double IM;

    ParticleVars(const TLorentzVector& lv, const ParticleTypeDatabase::Type& type) noexcept;
    ParticleVars(const data::Particle& p) noexcept;
    ParticleVars(double e=0.0, double theta=0.0, double phi=0.0, double im=0.0) noexcept:
        E(e), Theta(theta), Phi(phi), IM(im) {}
    ParticleVars(const ParticleVars&) = default;
    ParticleVars(ParticleVars&&) = default;
    ParticleVars& operator=(const ParticleVars&) =default;
    ParticleVars& operator=(ParticleVars&&) =default;
    virtual void SetBranches(TTree* tree, const std::string& prefix);
    virtual void Clear();
};

struct ParticleTools {

    /**
     * @brief Construct a string describing the particle dacay tree using the particle::PrintName()s
     * @param particles
     * @return
     */
    static std::string GetDecayString(const data::ParticleTree_t& particletree);

    static std::string GetDecayString(const ParticleTypeTree& particletypetree);

    /**
     * @brief SanitizeDecayString replaces all special characters by _
     * @param decaystring input decaystring
     * @return string suitable to be used as histogram name prefix
     */
    static std::string SanitizeDecayString(std::string decaystring);


    /** @brief Construct a string describing the production channel (e.g. gamma p -> p pi0 for pion production)
     * @param particle
     * @return
     */
    static std::string GetProductionChannelString(const data::ParticleTree_t& particletree);

    /**
     * @brief Find the first Particle of given type in particle list
     * @param type Type to search for
     * @param particles List to search in
     * @return The Particle found
     */
    static const data::ParticlePtr FindParticle(const ParticleTypeDatabase::Type& type, const data::ParticleList& particles);

    static const data::ParticlePtr FindParticle(const ParticleTypeDatabase::Type& type, const data::ParticleTree_t& particletree, size_t maxlevel);


    /**
     * @brief FillIMCombinations loops over all n tuples of given particles, builds sum and fills invariant mass
     * @param h histogram to be filled with Fill(invariant mass)
     * @param n multiplicity or number of particles drawn from particles
     * @param particles list of particles
     */
    static void FillIMCombinations(TH1* h, unsigned n, const data::ParticleList& particles);

    static bool SortParticleByName(const data::ParticlePtr& a, const data::ParticlePtr& b);

    static bool MatchByParticleName(const data::ParticlePtr& a, const ParticleTypeDatabase::Type& b);

};

}
}

}
