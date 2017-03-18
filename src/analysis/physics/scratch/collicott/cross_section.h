#pragma once

#include "physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "analysis/plot/HistogramFactory.h"
#include "analysis/plot/PromptRandomHist.h"

class TH1D;

namespace ant {
namespace analysis {
namespace utils {
/**
 * @brief The CrossSection class
 * This class stores a yield tree... stored as (particle, tagger channel) pair
 * Yield tree has:
 *  is MC = marker for MC
 *  data  = data string, could be decay from MC, or setup, user specified
 *  sp_X  = kinematics for scattered particle
 *  tc_X  = parameters for tagger channel
 *  missing_mass = missing mass (ignoring recoil)
 *
 * Note: prompt random fill weight is stored as tc_promptrandom
 * Yield->Draw("missing_mass","tc_promptrandom") produces pr-subtracted missing mass
 *
 * Also stores scaler tree:
 *  tagger_channel = dummy vector of tagger channel
 *  tagger_scalers = vector of scalers
 *  tagger_eff     = vector of tagging efficiency
 *  tagger_deff    = vector of tagging efficiency errors
 *  exp_livetime   = experimental livetime
 *
 * Scalers->Draw("tagger_eff:tagger_channel") produces 2d tagg eff vs. tagger channel plot
 * Scalers->Draw("tagger_eff[30]") produces 1d tagg eff values for channel 30
 *
 * USE:
 *
 * Once per event:
 *
 *      cross_section.SetEventType(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC), decay);
 *      Note: this will also track incident flux
 *
 * After all cuts,loop over tagger channels:
 *
 *       for (const auto &tc : event.Reconstructed().TaggerHits){
 *          promptrandom.SetTaggerHit(tc.Time);
 *          cross_section.AcceptEvent(Meson,Meson_time,tc,promptrandom);
 *       }
 *
 */

class scratch_collicott_CrossSection {
public:
    struct Yield_t : WrapTTree {
        Yield_t();

        // MC Decay branch
        ADD_BRANCH_T(bool, isMC)
        ADD_BRANCH_T(std::string, data)

        // SP kinematics
        ADD_BRANCH_T(double, sp_theta)
        ADD_BRANCH_T(double, sp_phi)
        ADD_BRANCH_T(double, sp_time)
        ADD_BRANCH_T(double, sp_im)

        // Tagger related kinematics
        ADD_BRANCH_T(double, tc_channel)
        ADD_BRANCH_T(double, tc_photonenergy)
        ADD_BRANCH_T(double, tc_time)
        ADD_BRANCH_T(double, tc_promptrandom)

        // Reaction related kinematics
        ADD_BRANCH_T(double, missing_mass)


    };

    struct Scalers_t : WrapTTree {
        Scalers_t();

        // Scaler branch
        ADD_BRANCH_T(std::vector<int>,      tagger_channel)
        ADD_BRANCH_T(std::vector<long>,     tagger_scalers)
        ADD_BRANCH_T(std::vector<double>,   tagger_eff)
        ADD_BRANCH_T(std::vector<double>,   tagger_deff)
        ADD_BRANCH_T(double,                exp_livetime)

    };

protected:

    std::shared_ptr<TaggerDetector_t> Tagger;
    unsigned int nTagger;

    HistogramFactory HistFac;


    TH1D*   flux;
    TH1D*   flux_ltcorrected;

    Yield_t yield;
    Scalers_t scalers;

    bool event_isMC = false;
    std::string event_decay = "";
    bool useSC = false;

    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
public:
    scratch_collicott_CrossSection(const HistogramFactory &histFac, OptionsPtr opts);
    virtual ~scratch_collicott_CrossSection() {}

    void SetEventType(bool isMC, const std::string& decay);
    void AcceptEvent(const LorentzVec &sp, double sp_time, const TTaggerHit &tc, const PromptRandom::Switch &promptrandom);
    void TrackIncidentFlux();

};

}
}
}
