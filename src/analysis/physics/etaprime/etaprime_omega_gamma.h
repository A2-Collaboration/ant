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

        ADD_BRANCH_T(std::vector<double>, ggg, 4)
        ADD_BRANCH_T(std::vector<double>, gg_gg1, 3)
        ADD_BRANCH_T(std::vector<double>, gg_gg2, 3)

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
        double          EBeam = std_ext::NaN;
    };

    struct Sig_t {

        // the subtree to be fitted is either pi0->2g
        // or omega->pi0g->3g
        // fitting the whole decay tree would overconstrain the
        // photons

        struct Fit_t {

            struct Tree_t : WrapTTree {

                ADD_BRANCH_T(double,   TreeFitChi2)
                ADD_BRANCH_T(double,   TreeFitProb)
                ADD_BRANCH_T(unsigned, TreeFitIterations)

                ADD_BRANCH_T(double, IM_Pi0_best)
                ADD_BRANCH_T(double, IM_Pi0_fitted)

                ADD_BRANCH_T(double, IM_Pi0gg)

                ADD_BRANCH_T(double, IM_gg)

                ADD_BRANCH_T(unsigned, MCTrueMatch)

                void Reset();
            };

            Fit_t(utils::TreeFitter fitter);

            static utils::TreeFitter Make(const ParticleTypeDatabase::Type& subtree);

            utils::TreeFitter treefitter;

            utils::TreeFitter::tree_t fitted_Pi0;
            utils::TreeFitter::tree_t fitted_g1_Pi0;
            utils::TreeFitter::tree_t fitted_g2_Pi0;

            utils::TreeFitter::tree_t fitted_Omega;
            utils::TreeFitter::tree_t fitted_g_Omega;

        };

        struct Pi0_t : Fit_t {

            Pi0_t();

            struct Tree_t : Fit_t::Tree_t {
                ADD_BRANCH_T(std::vector<double>, IM_Pi0g, 2)
                void Reset();
            };

            Tree_t t;

            void Process(const Particles_t& particles, TParticleTree_t ptree_sigref);
        };

        struct OmegaPi0_t : Fit_t {

            OmegaPi0_t();

            struct Tree_t : Fit_t::Tree_t {
                ADD_BRANCH_T(double, IM_Pi0g_fitted)
                ADD_BRANCH_T(double, IM_Pi0g_best)

                ADD_BRANCH_T(double, Bachelor_E_fitted)
                ADD_BRANCH_T(double, Bachelor_E_best)

                void Reset();
            };

            Tree_t t;

            void Process(const Particles_t& particles, TParticleTree_t ptree_sigref);

        };

        Pi0_t Pi0;
        OmegaPi0_t OmegaPi0;


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
