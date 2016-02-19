#pragma once

#include "analysis/physics/Physics.h"
#include "utils/particle_tools.h"
#include "utils/Fitter.h"
#include "plot/PromptRandomHist.h"

#include "base/ParticleTypeTree.h"

#include <cassert>

class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {
namespace physics {

class EtapOmegaG : public Physics {

    TH1D* h_CommonCuts;
    TH1D* h_MissedBkg;

    // variables for TTree branches
    // shared among sig/ref analyses

    TTree* treeCommon;
    bool     b_IsSignal;
    unsigned b_MCTrue;
    unsigned b_nPhotonsCB;
    unsigned b_nPhotonsTAPS;
    double   b_CBSumVetoE;
    double   b_CBAvgTime;
    double   b_PIDSumE;

    double   b_ProtonCopl;
    double   b_MissingMass;
    double   b_KinFitChi2;
    double   b_TaggW;
    double   b_TaggW_tight;
    double   b_TaggE;
    double   b_TaggT;
    unsigned b_TaggCh;

    PromptRandom::Switch promptrandom;
    PromptRandom::Switch promptrandom_tight;

    utils::KinFitter kinfitter_2;
    utils::KinFitter kinfitter_4;

    struct Particles_t {
        TParticlePtr    Proton;
        TParticleList   Photons;
        TLorentzVector  PhotonSum;
    };

    struct Sig_t {
        Sig_t();

        TTree* Tree;
        utils::TreeFitter treefitter;

        utils::TreeFitter::tree_t fitted_EtaPrime;
        utils::TreeFitter::tree_t fitted_Omega;
        utils::TreeFitter::tree_t fitted_Pi0;

        utils::TreeFitter::tree_t fitted_g_EtaPrime;
        utils::TreeFitter::tree_t fitted_g_Omega;
        utils::TreeFitter::tree_t fitted_g1_Pi0;
        utils::TreeFitter::tree_t fitted_g2_Pi0;


        double   b_TreeFitChi2;
        unsigned b_TreeFitIterations;

        double b_IM_EtaPrime_fitted;
        double b_IM_Omega_fitted;
        double b_IM_Pi0_fitted;

        double b_IM_EtaPrime_best;
        double b_IM_Omega_best;
        double b_IM_Pi0_best;


        double b_Bachelor_best_best;
        double b_Bachelor_best_fit;
        double b_Bachelor_fit_best;
        double b_Bachelor_fit_fit;

        unsigned b_MCTrueMatch;

        void SetupBranches();
        void ResetBranches();
        void Process(const Particles_t& particles, TParticleTree_t particletree);
    };

    struct Ref_t {
        TTree* Tree;
        double b_IM_2g;
        void SetupBranches();
        void ResetBranches();
        void Process(const Particles_t& particles);
    };

    Sig_t Sig;
    Sig_t SigFitted;
    Ref_t Ref;
    Ref_t RefFitted;

    friend class Sig_t;
    static const ParticleTypeTree ptreeSignal;
    static const ParticleTypeTree ptreeReference;
    static const std::vector<ParticleTypeTree> ptreeBackgrounds;

public:
    EtapOmegaG(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

class EtapOmegaG_MC : public Physics {

    struct expected_peak_t {
        double Mean;
        double Sigma;
        expected_peak_t(double mean, double sigma) :
            Mean(mean), Sigma(sigma) {}
        interval<double> makeCutInterval(unsigned nSigma=2) const {
            return {Mean-nSigma*Sigma,Mean+nSigma*Sigma};
        }
    };

    // means/sigma extracted from gg/ggg/gggg histograms for signal channel
    const expected_peak_t Pi0 = {126, 15};
    const expected_peak_t Eta = {515, 18};
    const expected_peak_t Omega = {735, 32};
    const expected_peak_t EtaPrime_sig = {895, 27};
    // extracted from gg histogram for reference channel
    const expected_peak_t EtaPrime_ref = {905, 29};

    ParticleTypeTree treeSignal;
    ParticleTypeTree treeReference;


    SmartHistFactory sig_HistFac;
    SmartHistFactory ref_HistFac;

    TH1D* h_TotalEvents;

    struct histogram_t {
        TH1D* Steps;
        TH1D* MissedBkg;
        histogram_t(SmartHistFactory HistFac);
    };

    histogram_t sig_hists;
    histogram_t ref_hists;

    template<typename T>
    struct perDecayHists_t {
        ParticleTypeTree Tree;
        std::string ShortName;
        std::string DecayString;
        T PerDecayHists;
        perDecayHists_t(const std::string& shortName,
                    const SmartHistFactory& HistFac_parent,
                    ParticleTypeTree tree
                    ) :
            Tree(tree),
            ShortName(shortName),
            DecayString(utils::ParticleTools::GetDecayString(tree)),
            PerDecayHists(SmartHistFactory(ShortName, HistFac_parent, DecayString))
        {}
        perDecayHists_t(const std::string& shortName, const SmartHistFactory& HistFac_parent) :
            Tree(nullptr),
            ShortName(shortName),
            DecayString(shortName),
            PerDecayHists(SmartHistFactory(ShortName, HistFac_parent, DecayString))
        {}
    };

    struct sig_perDecayHists_t {
        TH1D* Steps;

        TH1D* gggg;
        TH1D* ggg;
        TH1D* gg;

        TH2D* IM_gg_gg;
        TH2D* IM_gg_gg_cut;


        TH1D* Proton_Copl;

        TH2D* IM_etap_omega;
        TH1D* IM_pi0;

        TH1D* MM_gggg;
        TH1D* MM_etap;

        TH1D* Chi2_All;
        TH1D* Chi2_Best;

        TH1D* g_EtaPrime_E;

        sig_perDecayHists_t(SmartHistFactory HistFac);
    };

    std::vector<perDecayHists_t<sig_perDecayHists_t>> sig_perDecayHists;

    struct ref_perDecayHists_t {
        TH1D* Steps;

        TH1D* gg;

        TH2D* Proton_ThetaPhi;
        TH1D* Proton_Energy;

        TH1D* IM_etap;

        TH1D* MM_etap;

        TH1D* Proton_Copl;

        ref_perDecayHists_t(SmartHistFactory HistFac);
    };

    std::vector<perDecayHists_t<ref_perDecayHists_t>> ref_perDecayHists;

    struct sig_TTree_t {
        sig_TTree_t(TTree* tree) : Tree(tree) {}

        TTree* Tree;
        int MCTrueIndex = -1;

        utils::ParticleVars Proton;
        utils::ParticleVars ProtonTrue;

        double Chi2;
        utils::ParticleVars g_Pi0_0;
        utils::ParticleVars g_Pi0_1;
        utils::ParticleVars g_Omega;
        utils::ParticleVars g_EtaPrime;
        utils::ParticleVars g_EtaPrime_Boosted;
        utils::ParticleVars Pi0;
        utils::ParticleVars Omega;
        utils::ParticleVars EtaPrime;

        void SetBranches();
    };
    sig_TTree_t sig_TTree;

    struct ref_TTree_t {
        ref_TTree_t(TTree* tree) : Tree(tree) {}
        TTree* Tree;
        int MCTrueIndex = -1;

        utils::ParticleVars Proton;

        void SetBranches();
    };
    ref_TTree_t ref_TTree;


    void ProcessSig(const TParticleTree_t& particletree, const TEventData& data);
    void ProcessRef(const TParticleTree_t& particletree, const TEventData& data);


    template<typename T>
    const T& getHistogram(const TParticleTree_t& particletree,
                          const std::vector<EtapOmegaG_MC::perDecayHists_t<T>>& perDecayHists,
                          int& index
                          ) {
        assert(!perDecayHists.empty());
        index = perDecayHists.size()-1;
        if(!particletree)
            return perDecayHists.back().PerDecayHists;
        for(size_t i=0;i<perDecayHists.size()-1;i++) {
            auto& item = perDecayHists[i];
            if(particletree->IsEqual(item.Tree, utils::ParticleTools::MatchByParticleName)) {
                index = i;
                return item.PerDecayHists;
            }
        }
        return perDecayHists.back().PerDecayHists;
    }

public:
    EtapOmegaG_MC(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};


}}}
