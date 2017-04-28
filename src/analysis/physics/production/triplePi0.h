#pragma once

#include "analysis/utils/fitter/TreeFitter.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/TriggerSimulation.h"

#include "base/WrapTTree.h"

#include "TLorentzVector.h"

#include "analysis/physics/scratch/wolfes/tools/tools.h"


class TH1D;

namespace ant {
namespace analysis {
namespace physics {



struct triplePi0 :  Physics {

    //===================== Settings   ========================================================


    struct settings_t
    {
        enum class selectOn{
            kinFit,
            sigFit
        };
        selectOn selType = selectOn::kinFit;
        const std::string Tree_Name = "tree";

        const unsigned nPhotons = 6;

        const interval<size_t> Cut_NCands = {7,7};
        const interval<size_t> Cut_NNeutr = {5,7};
        const IntervalD Cut_ProtonCopl    = {-25,25};
        const IntervalD Cut_MM            = ParticleTypeDatabase::Proton.GetWindow(350).Round();
        const IntervalD Cut_MMAngle       = {0,25};
        const IntervalD Cut_EMB_prob      = {0.005,1};

        const IntervalD              Range_Prompt  =   { -5,  5};
        const std::vector<IntervalD> Ranges_Random = { {-55,-10},
                                                       { 10, 55}  };

        const double fitter_ZVertex = 3;
        const double vetoThreshE    = 0.0;

        const unsigned Index_Data     = 0;
        const unsigned Index_Signal   = 1;
        const unsigned Index_MainBkg  = 2;
        const unsigned Index_SigmaBkg = 3;
        const unsigned Index_Offset   = 10;
        const unsigned Index_Unknown  = 9;
    };

    static std::string getOtherChannelNames(const unsigned i);
    const settings_t phSettings;
    const std::shared_ptr<TaggerDetector_t> tagger;

    //===================== Channels   ========================================================

    struct named_channel_t
    {
        const std::string Name;
        const ParticleTypeTree DecayTree;
        named_channel_t(const std::string& name, ParticleTypeTree tree):
            Name(name),
            DecayTree(tree){}
    };

    static const named_channel_t              signal;
    static const named_channel_t              mainBackground;
    static const named_channel_t              sigmaBackground;
    static const std::vector<named_channel_t> otherBackgrounds;

    //===================== Histograms ========================================================

    TH1D* hist_steps             = nullptr;
    TH1D* hist_channels          = nullptr;
    TH1D* hist_channels_end      = nullptr;
    TH2D* hist_neutrals_channels = nullptr;


    //===================== KinFitting ========================================================


    std::shared_ptr<utils::UncertaintyModel> uncertModel;

    utils::KinFitter kinFitterEMB;

    utils::TreeFitter fitterSig;
    std::vector<utils::TreeFitter::tree_t> pionsFitterSig;

//    utils::TreeFitter fitterSigmaPlus;
//    std::vector<utils::TreeFitter::tree_t> pionsFitterSigmaPlus;
//    utils::TreeFitter::tree_t kaonFitterSigmaPlus;
//    utils::TreeFitter::tree_t sigmaFitterSigmaPlus;


    //========================  ProptRa. ============================================================

    utils::TriggerSimulation triggersimu;
    ant::analysis::PromptRandom::Switch promptrandom;

    //========================  Storage  ============================================================


    struct fitRatings_t
    {
        double Prob;
        double Chi2;
        int    Niter;
        bool   FitOk;
        TLorentzVector Proton;
        std::vector<TLorentzVector>   Intermediates;
        std::vector<unsigned>         PhotonCombination;
        fitRatings_t(double prob,double chi2,int niter, bool fitOk,
                     const TLorentzVector& proton,
                     const std::vector<TLorentzVector>&   intermediates,
                     const std::vector<unsigned>&         photonCombination):
            Prob(prob),Chi2(chi2),Niter(niter), FitOk(fitOk),
            Proton(proton),
            Intermediates(intermediates),PhotonCombination(photonCombination){}
    };

    struct PionProdTree : WrapTTree
    {

        // type: 0   data
        //       1   signal (3pi0)
        //       2   mainBkg(eta->3pi0)
        //       10+ otherBkg
        ADD_BRANCH_T(unsigned, MCTrue)

        ADD_BRANCH_T(double,   Tagg_W)
        ADD_BRANCH_T(unsigned, Tagg_Ch)
        ADD_BRANCH_T(double,   Tagg_E)
        ADD_BRANCH_T(double,   Tagg_Eff)
        ADD_BRANCH_T(double,   Tagg_EffErr)

        ADD_BRANCH_T(double,   ChargedClusterE)
        ADD_BRANCH_T(double,   ChargedCandidateE)
        ADD_BRANCH_T(unsigned, Neutrals)

        ADD_BRANCH_T(double, CBAvgTime)
        ADD_BRANCH_T(double, CBESum)

        ADD_BRANCH_T(double, PhotonVeto)
        ADD_BRANCH_T(double, ProtonVeto)
        ADD_BRANCH_T(std::vector<double>, CBPhotonVeto)

        // best emb combination raw
        ADD_BRANCH_T(TLorentzVector,              proton)
        ADD_BRANCH_T(double,                      protonTime)

        ADD_BRANCH_T(std::vector<TLorentzVector>, photons)
        ADD_BRANCH_T(std::vector<double>,         photonTimes)
        ADD_BRANCH_T(TLorentzVector,              photonSum)
        ADD_BRANCH_T(double,                      IM6g)

        ADD_BRANCH_T(TLorentzVector,              proton_MM)
        ADD_BRANCH_T(double,                      pMM_angle)
        ADD_BRANCH_T(double,                      pg_copl)
        void SetRaw(const tools::protonSelection_t& selection);

        // best emb comb. emb-fitted
        ADD_BRANCH_T(TLorentzVector,              EMB_proton)
        ADD_BRANCH_T(std::vector<TLorentzVector>, EMB_photons)
        ADD_BRANCH_T(TLorentzVector,              EMB_photonSum)
        ADD_BRANCH_T(double,                      EMB_IM6g)
        ADD_BRANCH_T(double,                      EMB_Ebeam)
        ADD_BRANCH_T(TLorentzVector,              EMB_proton_MM)

        ADD_BRANCH_T(double,                      EMB_prob)
        ADD_BRANCH_T(double,                      EMB_chi2)
        ADD_BRANCH_T(int,                         EMB_iterations)
        void SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result);

        //best tree-fit combination
        ADD_BRANCH_T(double,                      SIG_prob)
        ADD_BRANCH_T(double,                      SIG_chi2)
        ADD_BRANCH_T(int,                         SIG_iterations)
        ADD_BRANCH_T(std::vector<TLorentzVector>, SIG_pions)
        ADD_BRANCH_T(TLorentzVector,              SIG_proton)
        ADD_BRANCH_T(double,                      SIG_IM3Pi0)
        ADD_BRANCH_T(std::vector<unsigned>,       SIG_combination)
        ADD_BRANCH_T(double,                      SIG_photonVeto)
        ADD_BRANCH_T(double,                      SIG_corrPhotonVeto)
        void SetSIG(const triplePi0::fitRatings_t& fitRating);

        ADD_BRANCH_T(double,                      BKG_prob)
        ADD_BRANCH_T(double,                      BKG_chi2)
        ADD_BRANCH_T(int,                         BKG_iterations)
        ADD_BRANCH_T(std::vector<TLorentzVector>, BKG_pions)
        ADD_BRANCH_T(TLorentzVector,              BKG_proton)
        ADD_BRANCH_T(std::vector<unsigned>,       BKG_combination)
        void SetBKG(const triplePi0::fitRatings_t& fitRating);

        ADD_BRANCH_T(double,                      SIGMA_prob)
        ADD_BRANCH_T(double,                      SIGMA_chi2)
        ADD_BRANCH_T(int,                         SIGMA_iterations)
        ADD_BRANCH_T(std::vector<TLorentzVector>, SIGMA_pions)
        ADD_BRANCH_T(TLorentzVector,              SIGMA_k0s)
        ADD_BRANCH_T(TLorentzVector,              SIGMA_SigmaPlus)
        ADD_BRANCH_T(std::vector<unsigned>,       SIGMA_combination)

        static constexpr const char* treeName()       {return "tree";}
        static constexpr const char* treeAccessName() {return "triplePi0/tree";}
    };
    PionProdTree tree;

    //========================  MAIN     ============================================================

    triplePi0(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override {}
    virtual void ShowResult() override;


    void FillStep(const std::string& step) {hist_steps->Fill(step.c_str(),1);}

};


}}}
