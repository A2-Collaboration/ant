#pragma once

#include "tree/TParticle.h"
#include "tree/TTaggerHit.h"

#include <functional>

namespace ant {
namespace analysis {
namespace utils {

constexpr IntervalD nocut{-std_ext::inf, std_ext::inf};

/**
 * @brief The ProtonPhotonCombs struct provides convinient management of the proton/photons loop
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
        // set when filtered
        TTaggerHit TaggerHit{0, std_ext::NaN, std_ext::NaN};
        double MissingMass{std_ext::NaN};
        LorentzVec PhotonSum{{0,0,0},0};
        double DiscardedEk{std_ext::NaN};
    };

    using Observer_t =  std::function<void(const std::string&)>;

    /**
     * @brief The Combinations_t struct manages the available proton/photon combinations as a whole
     */
    struct Combinations_t : std::list<comb_t> {

        Combinations_t& FilterMult(unsigned nPhotonsRequired, double maxDiscardedEk = std_ext::inf);

        Combinations_t& FilterIM(const IntervalD& photon_IM_sum_cut = nocut);

        Combinations_t& FilterMM(const TTaggerHit& taggerhit, const IntervalD& missingmass_cut = nocut,
                                 const ParticleTypeDatabase::Type& target = ParticleTypeDatabase::Proton);
        Observer_t Observer;
    };

    // make them const which prevents accidental filtering and
    const Combinations_t Combinations;

    // use this operator to make a new copy, and observe filtering by passing lambda
    Combinations_t operator()(const Observer_t& observer = [] (const std::string&) {} ) const {
        auto copy = Combinations;
        copy.Observer = observer;
        return copy;
    }

    ProtonPhotonCombs(const TCandidateList& cands) :
        Combinations(MakeCombinations(cands))
    {} // empty ctor

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

protected:
    static Combinations_t MakeCombinations(const TCandidateList& cands);
};

}}} // namespace ant::analysis::utils