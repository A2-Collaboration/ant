#pragma once

#include "analysis/data/Particle.h"
#include <string>

class TH1;

namespace ant {
namespace analysis {
namespace utils {

struct ParticleTools {

    /**
     * @brief Construct a string describing the particle dacay tree using the particle::PrintName()s
     * @param particles
     * @return
     */
    static std::string GetDecayString(const data::ParticleTree_t& particletree);

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
    static const ant::analysis::data::ParticlePtr FindParticle(const ant::ParticleTypeDatabase::Type& type, const data::ParticleList& particles);


    /**
     * @brief FillIMCombinations loops over all n tuples of given particles, builds sum and fills invariant mass
     * @param h histogram to be filled with Fill(invariant mass)
     * @param n multiplicity or number of particles drawn from particles
     * @param particles list of particles
     */
    static void FillIMCombinations(TH1* h, unsigned n, const data::ParticleList& particles);

    static bool SortParticleByName(const data::ParticlePtr& a, const data::ParticlePtr& b);

    static bool SortParticleTypeByName(const ParticleTypeDatabase::Type& a, const ParticleTypeDatabase::Type& b);

    static bool MatchByParticleName(const data::ParticlePtr& a, const ParticleTypeDatabase::Type& b);

};

}
}

}
