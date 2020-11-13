#pragma once

#include "analysis/physics/Physics.h"
#include "utils/TriggerSimulation.h"
#include "plot/PromptRandomHist.h"
#include "utils/fitter/KinFitter.h"
#include "utils/fitter/NoProtonFitter.h"
#include "analysis/utils/Uncertainties.h"
#include "utils/Matcher.h"
#include "analysis/utils/ProtonPhotonCombs.h"
#include "expconfig/detectors/TAPS.h"
#include "expconfig/detectors/CB.h"
#include <string>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class scratch_lheijken_checkkinfit : public Physics {
protected:
    utils::TriggerSimulation triggersimu;
    PromptRandom::Switch promptrandom;
    std::shared_ptr<expconfig::detector::TAPS> taps_detector;
    std::shared_ptr<expconfig::detector::CB> cb_detector;
    using particle_comb_t = utils::ProtonPhotonCombs::comb_t;
    using particle_combs_t = utils::ProtonPhotonCombs::Combinations_t;
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    enum parttype{en_p=0, en_ep, en_em, en_g};
    static const int nrFitCases = 5;
    static const int nrPcut = 3;
    static const int nrPartTypes = 4;
    static const int nrKinVars = 3;
    static const int nrFitVars = 4;
    double TAPSZPos;
    double CBRad;
    static const int chooserecprot = 0; //0 - all events, 1 - only events with rec. proton, 2 - only events without rec. proton
    static const int chooseimreg = 0; //0 - all events, 1 - only events in "peak" region, 2 - only events in "bg" region

    utils::UncertaintyModelPtr fit_model;
    utils::KinFitter fitter;
    utils::NoProtonFitter npfitter;

    TH1D *h_Steps, *h_Probability[nrFitCases][nrPcut];
    TH2D *h_TrueRec_vsEk_CB[nrPartTypes][nrKinVars], *h_TrueRec_vsTh_CB[nrPartTypes][nrKinVars], *h_TrueRec_vsPh_CB[nrPartTypes][nrKinVars];
    TH2D *h_TrueRec_vsEk_TAPS[nrPartTypes][nrKinVars], *h_TrueRec_vsTh_TAPS[nrPartTypes][nrKinVars], *h_TrueRec_vsPh_TAPS[nrPartTypes][nrKinVars];
    TH2D *h_FitRec_vsEk_CB[nrFitCases][nrPcut][nrPartTypes][nrKinVars], *h_FitRec_vsTh_CB[nrFitCases][nrPcut][nrPartTypes][nrKinVars], *h_FitRec_vsPh_CB[nrFitCases][nrPcut][nrPartTypes][nrKinVars];
    TH2D *h_FitRec_vsEk_TAPS[nrFitCases][nrPcut][nrPartTypes][nrKinVars], *h_FitRec_vsTh_TAPS[nrFitCases][nrPcut][nrPartTypes][nrKinVars], *h_FitRec_vsPh_TAPS[nrFitCases][nrPcut][nrPartTypes][nrKinVars];
    TH2D *h_TrueFit_vsEk_CB[nrFitCases][nrPcut][nrPartTypes][nrKinVars], *h_TrueFit_vsTh_CB[nrFitCases][nrPcut][nrPartTypes][nrKinVars], *h_TrueFit_vsPh_CB[nrFitCases][nrPcut][nrPartTypes][nrKinVars];
    TH2D *h_TrueFit_vsEk_TAPS[nrFitCases][nrPcut][nrPartTypes][nrKinVars], *h_TrueFit_vsTh_TAPS[nrFitCases][nrPcut][nrPartTypes][nrKinVars], *h_TrueFit_vsPh_TAPS[nrFitCases][nrPcut][nrPartTypes][nrKinVars];
    TH2D *h_TrueFit_zvert[nrFitCases][nrPcut], *h_TrueRec_Ebeam, *h_FitRec_Ebeam[nrFitCases][nrPcut], *h_TrueFit_Ebeam[nrFitCases][nrPcut];
    TH1D *h_IMgg_True, *h_IMgg_Rec[nrFitCases][nrPcut], *h_IMgg_Fit[nrFitCases][nrPcut], *h_IMeeg_True, *h_IMeeg_Rec[nrFitCases][nrPcut], *h_IMeeg_Fit[nrFitCases][nrPcut];
    TH1D *h_PartPulls_CB[nrFitCases][nrPcut][nrPartTypes][nrFitVars], *h_PartPulls_TAPS[nrFitCases][nrPcut][nrPartTypes][nrFitVars], *h_EbeamPulls[nrFitCases][nrPcut], *h_ZvertPulls[nrFitCases][nrPcut];

    void CreateHistos();
    void DoMatchTrueRecoStuff(const TParticleList &allmcpart, const std::vector<TParticlePtr> &trueparts, const TCandidateList &recocands, std::vector<std::vector<TParticlePtr>> &matchtruerecpart, double zvert);
    void DoFitComparisons(const int kfnr, const int pcnr, const int nrphif, const LorentzVec fitprot, const TParticlePtr recprot, const std::vector<LorentzVec> fitphot, const TParticleList recphot, const std::vector<std::vector<TParticlePtr>> matchtruerecopart, const std::vector<TParticlePtr> truepart, const std::vector<utils::Fitter::FitParticle> fitparts, double zvertfit, double tw);

public:
    scratch_lheijken_checkkinfit(const std::string& name, OptionsPtr opts);

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned max_iterations);

    virtual void ProcessEvent(const TEvent& event, manager_t&) override;
    virtual void ShowResult() override;
};

}}}
