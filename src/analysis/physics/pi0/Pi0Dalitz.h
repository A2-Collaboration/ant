#pragma once

#include "analysis/physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "expconfig/detectors/Tagger.h"
#include "expconfig/detectors/PID.h"
#include "expconfig/detectors/TAPSVeto.h"
#include "expconfig/detectors/CB.h"
#include "analysis/utils/TriggerSimulation.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/ProtonPhotonCombs.h"
#include "analysis/utils/fitter/KinFitter.h"
#include "analysis/utils/fitter/NoProtonFitter.h"
#include "analysis/utils/Uncertainties.h"
#include "TLorentzVector.h"
#include "analysis/utils/ClusterECorr_simple.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A class for analyzing the pi0->e+e-g decay
 *
 */
class Pi0Dalitz: public Physics {

public:
    struct TFFTree_t : WrapTTree {
        ADD_BRANCH_T(bool, TBIsDalDec)
        ADD_BRANCH_T(bool, TBIs2gDec)
        ADD_BRANCH_T(bool, TBHasProtonCand)
        ADD_BRANCH_T(std::vector<double>, TBIMeegs)
        ADD_BRANCH_T(std::vector<double>, TBIMees)
        ADD_BRANCH_T(std::vector<double>, TBIMggs)
        ADD_BRANCH_T(std::vector<bool>, TBsamePIDs)
        ADD_BRANCH_T(std::vector<double>, TBKFprobs)
        ADD_BRANCH_T(std::vector<double>, TBTaggWeights)
    };


protected:
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    using particle_comb_t = utils::ProtonPhotonCombs::comb_t;
    using particle_combs_t = utils::ProtonPhotonCombs::Combinations_t;

    std::shared_ptr<expconfig::detector::Tagger> tagger_detector;
    std::shared_ptr<expconfig::detector::CB> cb_detector;
    std::shared_ptr<expconfig::detector::PID> pid_detector;
    std::shared_ptr<expconfig::detector::TAPSVeto> veto_detector;

    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    enum parttype{en_p=0, en_ep, en_em, en_g};
    static const int nrPartTypes = 4;
    static const bool doecorr = true;
    const double pidecut_chaneu = 0.2;
    const double vetoecut_chaneu = 0.;
    const double pidecut_lepcan = 3;
    const double pidcbcut_protcan_sl = -0.02;
    const double pidcbcut_protcan_in = 3;
    const double vetotacut_protcan_sl = -0.015;
    const double vetotacut_protcan_in = 3;

    utils::ClusterECorr_simple lepCBcorr;
    utils::ClusterECorr_simple photCBcorr;
    utils::ClusterECorr_simple photTAPScorr;

    utils::UncertaintyModelPtr fit_model;
    utils::KinFitter fitter;
    utils::NoProtonFitter npfitter;

    TFFTree_t TFFTree;

    void CreateHistos();
    void DoTrueMCStuff(const int cut, const std::vector<bool> &WhichMC, const std::vector<TParticlePtr> &trueparts, const double &tw);
    void DoCandSelStuff(const TCandidateList &recocands, TCandidatePtrList &selcands, std::vector<TParticlePtr> &selprot, std::vector<TParticleList> &selphots, std::vector<bool> &cases);
    void PlotCandSelStuff(const TCandidatePtr cand, const int sel);
    TCandidate ParticleECorr(const TCandidatePtr cand, const bool &photon, const bool &incb);
    void DoMatchTrueRecoStuff(const TParticleList &allmcpart, const std::vector<TParticlePtr> &trueparts, const TCandidatePtrList &recocands, std::vector<TParticlePtr> &matchrecopart);
    void DoTaggerStuff(const int cut, const TLorentzVector &g, const double &time, const double &cortime, const double &tw);
    void DoTriggerStuff(const int cut, const double &tw);
    void DoRecoCandStuff(const int cut, const TCandidatePtrList &recocands, const std::vector<TParticlePtr> &selprot, const std::vector<TParticleList> &selphots, const std::vector<TParticlePtr> &recmatparts, const std::vector<bool> &WhichMC, const TLorentzVector &ig, const double &tw);
    double DoKinFitStuffWProt(const int KFind, const std::vector<TParticlePtr> &selprot, const std::vector<TParticleList> &selphots, const TLorentzVector &ig, utils::KinFitter &fitobj, std::vector<TParticlePtr> &bestprobrec, std::vector<TParticlePtr> &bestprobfit, const double tw);
    double DoKinFitStuffNProt(const int KFind, const std::vector<TParticleList> &selphots, const TLorentzVector &ig, utils::NoProtonFitter &fitobj, TParticleList &fitphotons, LorentzVec &fitproton, const double tw);

    //-- Histograms
    static const int nrEvSel = 13;
    enum evsel{en_nocut=0, en_cbsum, en_tagge, en_1p3g, en_1p3g_pc, en_1p3gnp_pc, en_np3g, en_np3g_pc, en_1p2g, en_1p2g_pc, en_1p2gnp_pc, en_np2g, en_np2g_pc};
    const double evselhistoprobcut = 0.05;
    static const int nrCandSel = 12;
    static const int nrKF = 6;
    static const int nrPcut = 8;
    const double pcuts[nrPcut] = {0., 0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.5};
    static const int nrKFVars = 4;
    //--- MCtrue
    TH1D *h_IMeegTrue, *h_IMggTrue, *h_IMeeTrue, *h_RecoTrueMatch, *h_EktrueEkrec[nrPartTypes], *h_EktrueEkrec_gg, *h_ee_angle[nrEvSel], *h_gg_angle[nrEvSel];
    TH2D *h_RecoTrueAngle, *h_ThetavsEnergy_MCTrue[nrEvSel][nrPartTypes];
    //--- Trigger
    TH1D *h_CBEsum[nrEvSel];
    //--- Tagger
    TH1D *h_TaggTimeuw[nrEvSel], *h_TaggcorTimeuw[nrEvSel], *h_TaggTimeww[nrEvSel], *h_TaggcorTimeww[nrEvSel],*h_TaggPhEnww[nrEvSel], *h_TaggPhEnuw[nrEvSel];
    //--- Reconstructed candidate info during candidate selection
    TH2D *h_PIDEvsT_cs[nrCandSel], *h_TVetoEvsT_cs[nrCandSel], *h_CBEvsT_cs[nrCandSel], *h_TAPSEvsT_cs[nrCandSel], *h_CBEvsNrCr_cs[nrCandSel], *h_TAPSEvsNrCr_cs[nrCandSel];
    TH2D *h_TimeCBvsPID_cs[nrCandSel], *h_EnergyCBvsPID_cs[nrCandSel], *h_TimeTAPSvsTVeto_cs[nrCandSel], *h_EnergyTAPSvsTVeto_cs[nrCandSel];
    //--- Reconstructed candidate info during event selection
    TH2D *h_PIDMultUsed, *h_VetoMultUsed, *h_NrRecCand[nrEvSel];
    TH1D *h_IMeegReco[nrEvSel], *h_IMggReco[nrEvSel], *h_MMpReco[nrEvSel], *h_OpAngpphReco[nrEvSel];
    TH2D *h_PIDEvsT[nrEvSel], *h_TVetoEvsT[nrEvSel], *h_CBEvsT[nrEvSel], *h_TAPSEvsT[nrEvSel], *h_CBEvsNrCr[nrEvSel], *h_TAPSEvsNrCr[nrEvSel];
    TH2D *h_TimeCBvsPID[nrEvSel], *h_EnergyCBvsPID[nrEvSel], *h_TimeTAPSvsTVeto[nrEvSel], *h_EnergyTAPSvsTVeto[nrEvSel];
    TH2D *h_ThPhCBvsPID[nrEvSel], *h_ThPhTAPSvsTVeto[nrEvSel], *h_ThetavsEnergy[nrEvSel];
    TH2D *h_EnergyCBvsPID_RecMat[nrEvSel][nrPartTypes+1], *h_EnergyTAPSvsTVeto_RecMat[nrEvSel][nrPartTypes];
    TH2D *h_ThetavsEnergy_RecMat[nrEvSel][nrPartTypes], *h_CBEvsNrCr_RecMat[nrEvSel][nrPartTypes], *h_TAPSEvsNrCr_RecMat[nrEvSel][nrPartTypes];
    TH1D *h_OpAngpphReco_RecMat[nrEvSel];
    //--- Overview
    TH1D *h_AnalysisStat;
    TH2D *h_AnalysisStat_RecMat;
    //--- KinFit
    TH1D *h_Stat[nrKF], *h_Status[nrKF], *h_Chi2[nrKF], *h_NDF[nrKF], *h_Prob[nrKF][nrPcut], *h_Zv[nrKF][nrPcut], *h_IM3g[nrKF][nrPcut], *h_IM2g[nrKF][nrPcut];
    TH2D *h_EP[nrKF], *h_NIter[nrKF], *h_NFuncCall[nrKF];
    TH1D *h_PartPulls_CB[nrKF][nrPcut][nrPartTypes][nrKFVars], *h_PartPulls_TAPS[nrKF][nrPcut][nrPartTypes][nrKFVars], *h_EbeamPulls[nrKF][nrPcut], *h_ZvertPulls[nrKF][nrPcut];


public:
    Pi0Dalitz(const std::string& name, OptionsPtr opts);
    virtual ~Pi0Dalitz() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned max_iterations);

};

}
}
}
