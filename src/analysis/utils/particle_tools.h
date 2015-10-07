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
    static std::string GetDecayString(const data::ParticleList& particles);

    /**
     * @brief SanitizeDecayString replaces all special characters by _
     * @param decaystring input decaystring
     * @return string suitable to be used as histogram name prefix
     */
    static std::string SanitizeDecayString(std::string decaystring);


    /** @brief Construct a string describing the production channel (e.g. gamma p -> p pi0 for pion production)
     *        Searches for a Beam+Target pseudo partcile and uses it
     *        and it's daughters to construct the string.
     *        "???" is returned if no Beam+Target particle is found.
     *        This should normally be run on the Intermediates() particle list of MCTrue data.
     * @param particles
     * @return
     */
    static std::string GetProductionChannelString(const data::ParticleList& particles);

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
};

}
}
}
