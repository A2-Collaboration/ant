#include "det_eff.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/ParticleTypeTree.h"
#include "slowcontrol/SlowControlVariables.h"
#include "expconfig/ExpConfig.h"

#include "TCanvas.h"
#include "TH1D.h"

#include <iostream>
#include <memory>


/**
 * Cross section bits
 */

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::utils;

scratch_collicott_DetEff::Stats_t::Stats_t() {}
scratch_collicott_DetEff::DetEff_t::DetEff_t() {}

scratch_collicott_DetEff::scratch_collicott_DetEff(const HistogramFactory& histFac, OptionsPtr) :
    HistFac("DetectionEfficiency",histFac)
{
    stats.CreateBranches(HistFac.makeTTree("Stats"));
    total.CreateBranches(HistFac.makeTTree("Sig_Total"));
    accepted.CreateBranches(HistFac.makeTTree("Sig_Accepted"));
    contamination.CreateBranches(HistFac.makeTTree("Background"));
}

void scratch_collicott_DetEff::SetEventType(bool isSignal, const std::string decay)
{
    stats.isSignal = isSignal;
    stats.decay = decay;
    stats.Tree->Fill();
}

void scratch_collicott_DetEff::TrackSignalEvent(const LorentzVec& sp, const double& sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom)
{
    FillDetEff(total, sp, sp_time, tc, promptrandom);
}

void scratch_collicott_DetEff::AcceptSigEvent(const LorentzVec& sp, const double& sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom)
{
    FillDetEff(accepted, sp, sp_time, tc, promptrandom);
}

void scratch_collicott_DetEff::AcceptBkgEvent(const LorentzVec& sp, const double& sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom)
{
    FillDetEff(contamination, sp, sp_time, tc, promptrandom);
}

void scratch_collicott_DetEff::FillDetEff(DetEff_t& t, const LorentzVec& sp, const double& sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom)
{
    t.reaction = stats.decay;

    t.sp_theta = sp.Theta()*radtodeg;
    t.sp_phi   = sp.Phi()*radtodeg;
    t.sp_time  = sp_time;

    t.tc_channel      = tc.Channel;
    t.tc_photonenergy = tc.PhotonEnergy;
    t.tc_time         = tc.Time;
    t.tc_promptrandom = promptrandom.FillWeight();

    t.Tree->Fill();
}
