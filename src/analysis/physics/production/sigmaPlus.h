#pragma once

#include "analysis/utils/fitter/TreeFitter.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/TriggerSimulation.h"

#include "analysis/utils/uncertainties/FitterSergey.h"

#include "utils/ProtonPhotonCombs.h"

#include "base/WrapTTree.h"

#include "TLorentzVector.h"
#include "tree/TSimpleParticle.h"

#include "analysis/physics/scratch/wolfes/tools/tools.h"

#include "tree/TSimpleParticle.h"


class TH1D;

namespace ant {
namespace analysis {
namespace physics {



struct sigmaPlus :  Physics {

    //===================== Settings   ========================================================

    struct settings_t
    {
        const std::string Tree_Name = "tree";

        const unsigned nPhotons = 6;

        const interval<size_t> Cut_NCands = {7,10};
        const interval<size_t> Cut_NNeutr = {5,7};
        const IntervalD Cut_ProtonCopl    = {-25,25};
        const IntervalD Cut_MM            = ParticleTypeDatabase::Proton.GetWindow(350).Round();
        const IntervalD Cut_MMAngle       = {0,25};
        const IntervalD Cut_EMB_prob      = {0.01,1};

        const IntervalD              Range_Prompt  =   { -5,  5};
        const std::vector<IntervalD> Ranges_Random = { {-55,-10},
                                                       { 10, 55}  };

        const double fitter_ZVertex = 3;
        const double vetoThreshE    = 0.0;

        static constexpr unsigned Index_Data       = 0;
        static constexpr unsigned Index_Signal     = 1;
        static constexpr unsigned Index_MainBkg    = 2;
        static constexpr unsigned Index_SumMC      = 3;
        static constexpr unsigned Index_BkgMC      = 4;

        static constexpr unsigned Index_unregTree  = 8;
        static constexpr unsigned Index_brokenTree = 9;
        static constexpr unsigned Index_Offset     = 10;
    };

    static std::string getOtherChannelNames(const unsigned i);
    const settings_t phSettings;
    const bool flag_mc;
    const std::shared_ptr<TaggerDetector_t> tagger;

    LorentzVec cmsBoost(const double taggE, const LorentzVec& part) const
    {
        TLorentzVector bt(0,0, taggE, taggE + ParticleTypeDatabase::Proton.Mass());
        auto comsPart = part;
        comsPart.Boost(-bt.BoostVector());
        return comsPart;
    }

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
    struct fitRatings_t
    {

        double Prob;
        double Chi2;
        int    Niter;
        bool   FitOk;
        TSimpleParticle Proton;
        std::vector<TLorentzVector>   Intermediates;
        std::vector<unsigned>         PhotonCombination;
        fitRatings_t(double prob,double chi2,int niter, bool fitOk,
                     const TSimpleParticle& proton,
                     const std::vector<TLorentzVector>&   intermediates,
                     const std::vector<unsigned>&         photonCombination):
            Prob(prob),Chi2(chi2),Niter(niter), FitOk(fitOk),
            Proton(proton),
            Intermediates(intermediates),PhotonCombination(photonCombination){}
    };

    std::shared_ptr<utils::UncertaintyModel> uncertModelData = std::make_shared<utils::UncertaintyModels::FitterSergey>();
    std::shared_ptr<utils::UncertaintyModel> uncertModelMC = std::make_shared<utils::UncertaintyModels::FitterSergey>();

    utils::KinFitter fitterEMB;

    utils::TreeFitter fitter3Pi0;
    std::vector<utils::TreeFitter::tree_t> pionsFitter3Pi0;
    fitRatings_t apply3Pi0Fit(const utils::ProtonPhotonCombs::comb_t& protonSelection,
                              const double Ebeam);

    utils::TreeFitter fitterK0S;
    utils::TreeFitter::tree_t K0SK0S;
    std::list<utils::TreeFitter::tree_t> PionsK0S;
    std::vector<utils::TreeFitter::tree_t> pionsFitterK0S;
    std::pair<LorentzVec,fitRatings_t> applyK0SFit(const utils::ProtonPhotonCombs::comb_t& protonSelection,
                                                        const double Ebeam);

    //========================  ProptRa. ============================================================

    utils::TriggerSimulation triggersimu;
    ant::analysis::PromptRandom::Switch promptrandom;

    //========================  Storage  ============================================================

    struct effTree_t : WrapTTree
    {

        ADD_BRANCH_T(double, Egamma)
        ADD_BRANCH_T(int,    TaggerBin)
        ADD_BRANCH_T(double, CosThetaK0S)
        virtual ~effTree_t(){}
    };

    struct SeenTree : effTree_t
    {
        static constexpr const char* treeName()       {return "seen";}
        static constexpr const char* treeAccessName() {return "sigmaPlus/seen";}
    };
    SeenTree seenSignal;

    struct RecTree : effTree_t
    {
        static constexpr const char* treeName()       {return "rec";}
        static constexpr const char* treeAccessName() {return "sigmaPlus/rec";}
    };
    RecTree recSignal;

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

        ADD_BRANCH_T(double,   ExpLivetime)

        ADD_BRANCH_T(double,   ChargedClusterE)
        ADD_BRANCH_T(double,   ChargedCandidateE)
        ADD_BRANCH_T(unsigned, Neutrals)
        ADD_BRANCH_T(double,   ProtonVetoE)
        ADD_BRANCH_T(double,   PionPIDVetoE)

        ADD_BRANCH_T(double,   CBAvgTime)
        ADD_BRANCH_T(double,   CBESum)

        ADD_BRANCH_T(int,      NCands)

        // best combination raw
        ADD_BRANCH_T(TSimpleParticle,              proton)
        ADD_BRANCH_T(std::vector<TSimpleParticle>, photons)
        ADD_BRANCH_T(TLorentzVector,               photonSum)
        ADD_BRANCH_T(double,                       IM6g)
        ADD_BRANCH_T(double,                       proton_MM)
        ADD_BRANCH_T(double,                       DiscardedEk)
        void SetRaw(const utils::ProtonPhotonCombs::comb_t& selection);

        // best comb. emb-fitted
        ADD_BRANCH_T(TSimpleParticle,              EMB_proton)
        ADD_BRANCH_T(std::vector<TSimpleParticle>, EMB_photons)
        ADD_BRANCH_T(TLorentzVector,               EMB_photonSum)
        ADD_BRANCH_T(double,                       EMB_IM6g)
        ADD_BRANCH_T(double,                       EMB_Ebeam)
        ADD_BRANCH_T(double,                       EMB_prob)
        ADD_BRANCH_T(double,                       EMB_chi2)
        ADD_BRANCH_T(int,                          EMB_iterations)
        void SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result);

        //best tree-fit combination
        ADD_BRANCH_T(TSimpleParticle,              SIG_proton)
        ADD_BRANCH_T(std::vector<TSimpleParticle>, SIG_photons)
        ADD_BRANCH_T(TLorentzVector,               SIG_photonSum)
        ADD_BRANCH_T(std::vector<TLorentzVector>,  SIG_pions)
        ADD_BRANCH_T(std::vector<unsigned>,        SIG_combination)
        ADD_BRANCH_T(double,                       SIG_IM6g)
        ADD_BRANCH_T(double,                       SIG_prob)
        ADD_BRANCH_T(double,                       SIG_chi2)
        ADD_BRANCH_T(int,                          SIG_iterations)
        void SetSIG(const sigmaPlus::fitRatings_t& fitRating);

        ADD_BRANCH_T(double,                       K0S_prob)

        ADD_BRANCH_T(TLorentzVector,               K0S)
        ADD_BRANCH_T(double,                       K0S_cosTheta)
        ADD_BRANCH_T(TLorentzVector,               K0S_Sigma)

        ADD_BRANCH_T(TLorentzVector,               K0S_EMB)
        ADD_BRANCH_T(double,                       K0S_cosTheta_EMB)
        ADD_BRANCH_T(TLorentzVector,               K0S_Sigma_EMB)

        ADD_BRANCH_T(double,                       K0S_IM6g)
        ADD_BRANCH_T(TLorentzVector,               K0S_photonSum)
        ADD_BRANCH_T(std::vector<TLorentzVector>,  K0S_pions)


        static constexpr const char* treeName()       {return "tree";}
        static constexpr const char* treeAccessName() {return "sigmaPlus/tree";}
    };
    PionProdTree tree;

    //========================  MAIN     ============================================================

    sigmaPlus(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override {}
    virtual void ShowResult() override;

    void FillStep(const std::string& step) {hist_steps->Fill(step.c_str(),1);}
};

}}}
