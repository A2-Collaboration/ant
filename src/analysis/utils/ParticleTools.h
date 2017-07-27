#pragma once

#include "Combinatorics.h"

#include "tree/TParticle.h"
#include "base/ParticleTypeTree.h"
#include "base/vec/LorentzVec.h"
#include "base/std_ext/misc.h"

#include "TH1.h"

#include <string>

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

    ParticleVars(const LorentzVec& lv, const ParticleTypeDatabase::Type& type) noexcept;
    ParticleVars(const TParticle& p) noexcept;
    ParticleVars(double e=0.0, double theta=0.0, double phi=0.0, double im=0.0) noexcept:
        E(e), Theta(theta), Phi(phi), IM(im) {}
    ParticleVars(const ParticleVars&) = default;
    ParticleVars(ParticleVars&&) = default;
    ParticleVars& operator=(const ParticleVars&) =default;
    ParticleVars& operator=(ParticleVars&&) =default;
    void SetBranches(TTree* tree, const std::string& prefix);
    void Clear();
};

class ParticleID;

struct ParticleTypeList {

    // factory to get a particle list from event using user-chosen particle ID
    static ParticleTypeList Make(const TCandidateList& cands);
    static ParticleTypeList Make(const TCandidateList& cands, const ParticleID& id);
    static ParticleTypeList Make(const TParticleTree_t& tree);
    static ParticleTypeList Make(const TParticleList& particles);

          TParticleList  GetAll() const && { return GetAll(); }
    const TParticleList& GetAll() const &  { return all; }

          TParticleList  Get(const ant::ParticleTypeDatabase::Type& type) const && { return Get(type); }
    const TParticleList& Get(const ant::ParticleTypeDatabase::Type& type) const & {
        auto entry = lists.find(std::addressof(type));
        if(entry == lists.end()) {
            static const TParticleList empty;
            return empty;
        }
        return entry->second;
    }

private:
    // force using the factory
    ParticleTypeList() = default;

    void Add(const TParticlePtr& particle)
    {
        lists[std::addressof(particle->Type())].emplace_back(particle);
        all.emplace_back(particle);
    }

    TParticleList all;
    std::map<const ParticleTypeDatabase::Type*, TParticleList> lists;
}; // ParticleTypeList

struct ParticleTools {

    /**
     * @brief GetPlutoProduction builds Pluto-compatible string of produced particles
     * @param particletypetree
     * @return representation suitable for Pluto's PReaction
     */
    static std::string GetPlutoProduction( const ParticleTypeTree& particletypetree);


    /**
     * @brief Construct a string describing the particle dacay tree using the particle::PrintName()s
     * @param particles
     * @return
     */
    static std::string GetDecayString(const TParticleTree_t& particletree, bool usePrintName = true);

    static std::string GetDecayString(const ParticleTypeTree& particletypetree, bool usePrintName = true);

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
    static std::string GetProductionChannelString(const TParticleTree_t& particletree);

    static ParticleTypeTree GetProducedParticle(const ParticleTypeTree& particletypetree);


    /**
     * @brief Find the first Particle of given type in particle list
     * @param type Type to search for
     * @param particles List to search in
     * @return The Particle found
     */
    static const TParticlePtr FindParticle(const ParticleTypeDatabase::Type& type, const TParticleList& particles);

    static const TParticlePtr FindParticle(const ParticleTypeDatabase::Type& type, const TParticleTree_t& particletree,
                                           size_t maxlevel = std::numeric_limits<size_t>::max());

    static const TParticleList FindParticles(const ParticleTypeDatabase::Type& type, const TParticleTree_t& particletree,
                                             size_t maxlevel = std::numeric_limits<size_t>::max());


    /**
     * @brief FillIMCombinations loops over all n tuples of given particle types, builds sum and fills invariant mass
     * @param h histogram to be filled with Fill(invariant mass)
     * @param n multiplicity or number of particles drawn from particles
     * @param particles list of particle types
     */
    template<typename T>
    static void FillIMCombinations(TH1* h, unsigned n, const std::vector<T>& particles, double weight = 1.0)
    {
        FillIMCombinations([h,weight] (double x) {h->Fill(x, weight);}, n, particles);
    }

    template<typename T>
    static void FillIMCombinations(std::function<void(double)> filler, unsigned n, const std::vector<T>& particles)
    {
        for(auto comb = makeCombination(particles,n); !comb.done(); ++comb) {
             LorentzVec sum({0,0,0},0);
             for(const auto& p : comb) {
                 sum += std_ext::dereference(p);
             }
             filler(sum.M());
        }
    }

    template<typename T>
    static void FillIMCombinations(std::vector<double>::iterator it, unsigned n, const std::vector<T>& particles)
    {
        FillIMCombinations([&it] (double v) {
            *it++ = v;
        }, n, particles);
    }


    static bool SortParticleByName(const TParticlePtr& a, const TParticlePtr& b);

    static bool MatchByParticleName(const TParticlePtr& a, const ParticleTypeDatabase::Type& b);

    static bool TryFindParticleDatabaseChannel(
            const TParticleTree_t& ptree,
            ParticleTypeTreeDatabase::Channel& channel);

};

}
}

}
