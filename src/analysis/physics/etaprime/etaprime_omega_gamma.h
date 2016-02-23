#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/utils/Fitter.h"
#include "analysis/plot/PromptRandomHist.h"

#include "base/ParticleTypeTree.h"
#include "base/WrapTTree.h"

#include <cassert>

class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {
namespace physics {

struct EtapOmegaG : Physics {

    TH1D* h_CommonCuts;
    TH1D* h_MissedBkg;

    // TreeCommon contains things
    // shared among sig/ref analyses

    struct TreeCommon : WrapTTree {
        ADD_BRANCH_T(bool,     IsSignal)
        ADD_BRANCH_T(unsigned, MCTrue)
        ADD_BRANCH_T(unsigned, nPhotonsCB)
        ADD_BRANCH_T(unsigned, nPhotonsTAPS)
        ADD_BRANCH_T(double,   CBSumE)
        ADD_BRANCH_T(double,   CBSumVetoE)
        ADD_BRANCH_T(double,   CBAvgTime)
        ADD_BRANCH_T(double,   ProtonTime)
        ADD_BRANCH_T(double,   PIDSumE)

        ADD_BRANCH_T(double,   ProtonCopl)
        ADD_BRANCH_T(double,   MissingMass)
        ADD_BRANCH_T(double,   KinFitChi2)
        ADD_BRANCH_T(unsigned, KinFitIterations)
        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggW_tight)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)
    };

    TreeCommon t;

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

        struct Tree_t : WrapTTree {

            ADD_BRANCH_T(std::vector<double>, ggg)
            ADD_BRANCH_T(std::vector<double>, gg_gg1)
            ADD_BRANCH_T(std::vector<double>, gg_gg2)

            ADD_BRANCH_T(double,   TreeFitChi2)
            ADD_BRANCH_T(unsigned, TreeFitIterations)

            ADD_BRANCH_T(double, IM_EtaPrime_fitted)
            ADD_BRANCH_T(double, IM_Omega_fitted)
            ADD_BRANCH_T(double, IM_Pi0_fitted)

            ADD_BRANCH_T(double, IM_EtaPrime_best)
            ADD_BRANCH_T(double, IM_Omega_best)
            ADD_BRANCH_T(double, IM_Pi0_best)

            ADD_BRANCH_T(double, Bachelor_best_best)
            ADD_BRANCH_T(double, Bachelor_best_fit)
            ADD_BRANCH_T(double, Bachelor_fit_best)
            ADD_BRANCH_T(double, Bachelor_fit_fit)

            ADD_BRANCH_T(unsigned, MCTrueMatch)

        };

        // we have multiple ideas for treefitting...
        struct Fit_t {

            Fit_t(const ParticleTypeDatabase::Type* typeptr = nullptr);

            static utils::TreeFitter MakeFitter(const ParticleTypeDatabase::Type* typeptr);

            utils::TreeFitter treefitter;

            utils::TreeFitter::tree_t fitted_EtaPrime;
            utils::TreeFitter::tree_t fitted_Omega;
            utils::TreeFitter::tree_t fitted_Pi0;

            utils::TreeFitter::tree_t fitted_g_EtaPrime;
            utils::TreeFitter::tree_t fitted_g_Omega;
            utils::TreeFitter::tree_t fitted_g1_Pi0;
            utils::TreeFitter::tree_t fitted_g2_Pi0;

            Tree_t t;

            void ResetBranches();
            void Process(const Particles_t& particles, TParticleTree_t particletree);

        };

        Sig_t();

        Fit_t All;
        Fit_t No_Pi0;
        Fit_t No_Omega;
        Fit_t No_EtaPrime;

        void SetupTrees(HistogramFactory HistFac);
        void Fill();
        void ResetBranches();
        void Process(const Particles_t& particles, TParticleTree_t particletree);


    };

    struct Ref_t {
        struct Tree_t : WrapTTree {
            ADD_BRANCH_T(double, IM_2g)
        };
        Tree_t t;

        void ResetBranches();
        void Process(const Particles_t& particles);
    };

    Sig_t Sig;
    Sig_t SigKinFit;
    Ref_t Ref;
    Ref_t RefKinFit;

    static const ParticleTypeTree ptreeSignal;
    static const ParticleTypeTree ptreeReference;
    struct Background_t {
        const std::string Name;
        const ParticleTypeTree Tree;
        Background_t(const std::string& name, ParticleTypeTree tree) :
            Name(name), Tree(tree) {}
    };
    static const std::vector<Background_t> ptreeBackgrounds;

    EtapOmegaG(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
