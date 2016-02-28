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
        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggW_tight)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)
    };

    TreeCommon t;

    PromptRandom::Switch promptrandom;
    PromptRandom::Switch promptrandom_tight;

    struct Particles_t {
        TParticlePtr    Proton;
        TParticleList   Photons;
        TLorentzVector  PhotonSum;
        double          EBeam = std_ext::NaN;
    };

    struct Sig_t {

        struct Tree_t : WrapTTree {

            ADD_BRANCH_T(std::vector<double>, ggg, 4)
            ADD_BRANCH_T(std::vector<double>, gg_gg1, 3)
            ADD_BRANCH_T(std::vector<double>, gg_gg2, 3)

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
            struct IM_Sigma_t {
                double EtaPrime;
                double Omega;
                IM_Sigma_t(double etaPrime = 1.0, double omega = 1.0) :
                    EtaPrime(etaPrime), Omega(omega) {}
                // Pi0 width is "reference"
            };

            Fit_t(utils::TreeFitter fitter);

            static utils::TreeFitter Make(const Fit_t::IM_Sigma_t& IM_Sigma,
                                          const ParticleTypeDatabase::Type* typeptr = nullptr);

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
            void Process(const Particles_t& particles, TParticleTree_t ptree_sigref);
            void DoPhotonCombinatorics(const TParticleList& photons);
            void CheckMCPhotonAssignment(const TParticleList& photons,
                                         TParticleTree_t ptree_sigref,
                                         TParticlePtr g_Omega_best,
                                         TParticlePtr g_EtaPrime_best);
        };

        struct FitOmegaPi0_t : Fit_t {
            FitOmegaPi0_t(utils::TreeFitter fitter);
            static utils::TreeFitter Make(const Fit_t::IM_Sigma_t& IM_Sigma);
            void Process(const Particles_t& particles, TParticleTree_t ptree_sigref);
        };



        Sig_t(const Fit_t::IM_Sigma_t& IM_Sigma = {});

        Fit_t All;
        Fit_t No_EtaPrime;
        FitOmegaPi0_t OmegaPi0;

        void SetupTrees(HistogramFactory HistFac);
        void Fill();
        void ResetBranches();
        void Process(const Particles_t& particles, TParticleTree_t ptree_sigref);


    };

    struct Ref_t {

        Ref_t();

        utils::KinFitter kinfitter;

        struct Tree_t : WrapTTree {
            ADD_BRANCH_T(double,   IM_2g)
            ADD_BRANCH_T(double,   KinFitChi2)
            ADD_BRANCH_T(unsigned, KinFitIterations)
        };
        Tree_t t;

        void ResetBranches();
        void Process(const Particles_t& particles);
    };

    Sig_t Sig;
    Ref_t Ref;

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
