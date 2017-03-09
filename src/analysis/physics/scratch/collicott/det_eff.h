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
 * @brief The DetEff class
 *
 * This class stores a set of yield trees needed to calculate detection efficiencies
 *
 * DetEff tree used to store full kinematics of
 * - total signal events (users should fill trueMC)
 * - accepted signal events (users should fill reconstructed)
 * - Store accepted background events (users should fill reconstructed)
 * - Basic stats
 *
 * USE:
 *
 * Once per event:
 *
 *       if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC))
 *               detection_efficiency.SetEventType(signal,decay);
 *
 * If signal event, and MC event, loop over tagger channels:
 *
 *       for (const auto &tc : event.Reconstructed().TaggerHits){
 *          promptrandom.SetTaggerHit(tc.Time);
 *          detection_efficiency.TrackSignalEvent(*pi0, 0, tc, promptrandom);
 *       }
 *
 * After all cuts, if MC event, loop over tagger channels:
 *
 *       for (const auto &tc : event.Reconstructed().TaggerHits){
 *          promptrandom.SetTaggerHit(tc.Time);
 *          detection_efficiency.AcceptEvent(Meson,Meson_time, tc,promptrandom);
 *       }
 *
 *
 */
class scratch_collicott_DetEff {
protected:

    struct Stats_t : WrapTTree {
        Stats_t();

        // SP kinematics
        ADD_BRANCH_T(bool, isSignal)
        ADD_BRANCH_T(std::string, decay)

    };

    struct DetEff_t : WrapTTree {
        DetEff_t();

        ADD_BRANCH_T(std::string, reaction)

        // SP kinematics
        ADD_BRANCH_T(double, sp_theta)
        ADD_BRANCH_T(double, sp_phi)
        ADD_BRANCH_T(double, sp_time)

        // Tagger related kinematics
        ADD_BRANCH_T(double, tc_channel)
        ADD_BRANCH_T(double, tc_photonenergy)
        ADD_BRANCH_T(double, tc_time)
        ADD_BRANCH_T(double, tc_promptrandom)
    };

    HistogramFactory HistFac;
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

    Stats_t  stats;
    DetEff_t total;
    DetEff_t accepted;
    DetEff_t contamination;

    bool Signal = false;
    std::string Decay = "";


public:
    scratch_collicott_DetEff(const HistogramFactory &histFac, OptionsPtr opts);
    virtual ~scratch_collicott_DetEff() {}

    void SetEventType(bool isSignal, const std::string decay);

    void TrackSignalEvent(const LorentzVec& sp, const double& sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom);
    void AcceptEvent(const LorentzVec& sp, const double& sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom)
    {
        if (stats.isSignal) AcceptSigEvent(sp,sp_time,tc,promptrandom);
        else                AcceptBkgEvent(sp,sp_time,tc,promptrandom);
    }

    void AcceptSigEvent( const LorentzVec& sp, const double& sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom);
    void AcceptBkgEvent(const LorentzVec& sp, const double& sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom);
    void FillDetEff(DetEff_t &t, const LorentzVec& sp, const double& sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom);

};

}
}
}
