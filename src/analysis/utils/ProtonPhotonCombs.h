#pragma once

#include "tree/TParticle.h"
#include "tree/TTaggerHit.h"

#include <functional>

namespace ant {
namespace analysis {
namespace utils {

// helper for defining default cut (which is no cut at all actually)
constexpr IntervalD nocut{-std_ext::inf, std_ext::inf};

/**
 * @brief The ProtonPhotonCombs struct provides convenient management of the proton/photons fitter loop
 */
struct ProtonPhotonCombs {
    /**
     * @brief The comb_t struct holds the information of one particular proton/photons selection
     *
     * The fields TaggerHit/MissingMass/PhotonSum/DiscardedEk are set by the Combinations::Filter* methods
     */
    struct comb_t {
        explicit comb_t(const TParticlePtr& proton) :
            Proton(proton) {}
        TParticleList Photons;
        TParticlePtr  Proton;

        double DiscardedEk{std_ext::NaN}; /// set by Combinations_t::FilterMult
        LorentzVec PhotonSum{{0,0,0},0};  /// set by Combinations_t::FilterIM
        double MissingMass{std_ext::NaN}; /// set by Combinations_t::FilterMM
    };

    using Observer_t =  std::function<void(const std::string&)>;

    /**
     * @brief The Combinations_t struct manages the available proton/photon combinations as a whole
     */
    struct Combinations_t : std::list<comb_t> {

        /**
         * @brief Observe sets the filtering observer and an optional prefix,
         * usually used to fill "steps" histogram
         * @param observer the observer, called for passes of filter
         * @param prefix always prepended to Observer() calls, keep it short
         * @return reference to modified instance
         */
        Combinations_t& Observe(const Observer_t& observer, const std::string& prefix = "") noexcept;

        /**
         * @brief FilterMult ensure multiplicity of photons is nPhotonsRequired,
         * rest are discarded (lowest Ek photons first, summed up in DiscardedEk)
         * @param nPhotonsRequired
         * @param maxDiscardedEk cut on DiscardedEk
         * @return reference to modified instance
         */
        Combinations_t& FilterMult(unsigned nPhotonsRequired, double maxDiscardedEk = std_ext::inf) noexcept;

        /**
         * @brief FilterIM kicks out combinations based on sum of photons invariant mass
         * @param photon_IM_sum_cut the allowed invariant mass (use std_ext::inf to define right-open interval)
         * @return reference to modified instance
         * @note should be called before FilterMM
         */
        Combinations_t& FilterIM(const IntervalD& photon_IM_sum_cut = nocut) noexcept;

        /**
         * @brief FilterMM kicks out proton/photons combinations based on missing mass
         * @param taggerhit the taggerhit to calculate the photon beam from
         * @param missingmass_cut the missing mass of the photons (should be around target type rest mass)
         * @param target the assumed target type (at rest), defaults to proton
         * @return reference to modified instance
         * @note FilterIM is automatically called if needed
         */
        Combinations_t& FilterMM(const TTaggerHit& taggerhit, const IntervalD& missingmass_cut = nocut,
                                 const ParticleTypeDatabase::Type& target = ParticleTypeDatabase::Proton) noexcept;

        using cut_t = std::function<bool(const comb_t&)>;
        /**
         * @brief FilterCustom specify own filter cut
         * @param cut function returning true if comb_t should be kicked out
         * @param name optional name for observing
         * @return reference to modified instance
         */
        Combinations_t& FilterCustom(const cut_t& cut,
                                     const std::string& name = "");

    private:
        Observer_t  Observer;
        std::string ObserverPrefix;
        bool called_FilterIM = false;
    };


    /**
     * @brief operator() call this to get copy of combinations for filtering (see above)
     * @return copy of pre-build combinations
     */
    Combinations_t operator()() const noexcept { return Combinations; }

    /**
     * @brief ProtonPhotonCombs pre-builds the particle combinations from given candidates
     * @param cands typically pass event.Reconstructed().Candidates
     * @note call only once per ProcessEvent to stay performant
     */
    ProtonPhotonCombs(const TCandidateList& cands) :
        Combinations(MakeCombinations(cands))
    {} // empty ctor

private:
    const Combinations_t Combinations;
    static Combinations_t MakeCombinations(const TCandidateList& cands) noexcept;
};

}}} // namespace ant::analysis::utils