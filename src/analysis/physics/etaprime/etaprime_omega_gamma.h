#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/utils/Fitter.h"
#include "analysis/utils/MCSmear.h"
#include "analysis/plot/PromptRandomHist.h"

#include "base/ParticleTypeTree.h"
#include "base/WrapTTree.h"

#include <cassert>

class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {

namespace utils {
class MCFakeReconstructed;
}

namespace physics {

struct EtapOmegaG : Physics {

    TH1D* h_CommonCuts;
    TH1D* h_CommonCuts_sig;
    TH1D* h_CommonCuts_ref;

    TH1D* h_MissedBkg;

    // TreeCommon contains things
    // shared among sig/ref analyses

    struct TreeCommon : WrapTTree {
        ADD_BRANCH_T(unsigned, MCTrue)
        ADD_BRANCH_T(double,   CBSumE)
        ADD_BRANCH_T(double,   CBAvgTime)
        ADD_BRANCH_T(double,   PIDSumE)

        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggW_tight)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)
    };

    struct SharedTree_t : WrapTTree {
        ADD_BRANCH_T(unsigned, nCandidates)
        ADD_BRANCH_T(double,   PhotonsEk)
        ADD_BRANCH_T(double,   DiscardedEk)

        ADD_BRANCH_T(unsigned, nPhotonsCB)
        ADD_BRANCH_T(unsigned, nPhotonsTAPS)
        ADD_BRANCH_T(double,   CBSumVetoE)

        ADD_BRANCH_T(double,   ProtonTime)
        ADD_BRANCH_T(double,   ProtonE)
        ADD_BRANCH_T(double,   ProtonTheta)
        ADD_BRANCH_T(double,   ProtonVetoE)
        ADD_BRANCH_T(double,   ProtonShortE)
        ADD_BRANCH_T(double,   ProtonTrueAngle)

        ADD_BRANCH_T(double,   PhotonSum)
        ADD_BRANCH_T(double,   ProtonCopl)
        ADD_BRANCH_T(double,   MissingMass)
        ADD_BRANCH_T(double,   FittedProtonE)

        ADD_BRANCH_T(double,   KinFitProb)
        ADD_BRANCH_T(int,      KinFitIterations)

        ADD_BRANCH_T(double,   KinFitBeamEPull)
        ADD_BRANCH_T(double,   KinFitProtonEPull)
        ADD_BRANCH_T(double,   KinFitProtonThetaPull)
        ADD_BRANCH_T(double,   KinFitProtonPhiPull)
        ADD_BRANCH_T(std::vector<double>,  KinFitPhotonEPulls)
        ADD_BRANCH_T(std::vector<double>,  KinFitPhotonThetaPulls)
        ADD_BRANCH_T(std::vector<double>,  KinFitPhotonPhiPulls)
    };

    TreeCommon t;

    PromptRandom::Switch promptrandom;
    PromptRandom::Switch promptrandom_tight;

    static const utils::UncertaintyModelPtr fit_uncertainty_model;

    utils::KinFitter kinfitter_sig;
    utils::KinFitter kinfitter_ref;

    std::unique_ptr<utils::MCSmear>             mc_smear;
    std::unique_ptr<utils::MCFakeReconstructed> mc_fake;


    struct Particles_t {
        double         PhotonEnergy;
        TParticlePtr   Proton;
        TParticleList  Photons;
        TParticleList  FittedPhotons;
        LorentzVec     PhotonSum;
        LorentzVec     FittedPhotonSum;
    };

    static bool findParticles(const TCandidatePtrList& candidates,
                              unsigned nPhotons,
                              TParticlePtr true_proton,
                              Particles_t& particles,
                              SharedTree_t& t,
                              TH1D* h_CommonCuts
                              );

    static bool doKinfit(const TTaggerHit& taggerhit,
                         utils::KinFitter& kinfitter,
                         Particles_t& particles,
                         SharedTree_t& t
                         );

    struct Sig_t {

        // the subtree to be fitted is either pi0->2g
        // or omega->pi0g->3g
        // fitting the whole decay tree would overconstrain the
        // photons

        struct Fit_t {

            struct BaseTree_t : WrapTTree {

                ADD_BRANCH_T(double, TreeFitProb)
                ADD_BRANCH_T(int,    TreeFitIterations)

                ADD_BRANCH_T(double,   TreeFitBeamEPull)
                ADD_BRANCH_T(double,   TreeFitProtonEPull)
                ADD_BRANCH_T(double,   TreeFitProtonThetaPull)
                ADD_BRANCH_T(double,   TreeFitProtonPhiPull)
                ADD_BRANCH_T(std::vector<double>,  TreeFitPhotonEPulls)
                ADD_BRANCH_T(std::vector<double>,  TreeFitPhotonThetaPulls)
                ADD_BRANCH_T(std::vector<double>,  TreeFitPhotonPhiPulls)

                ADD_BRANCH_T(double, IM_Pi0_best)
                ADD_BRANCH_T(double, IM_Pi0_fitted)

                ADD_BRANCH_T(double, IM_Pi0gg)

                ADD_BRANCH_T(double, IM_gg)

                ADD_BRANCH_T(unsigned, MCTrueMatch)

                void Reset();
            };

        protected:

            Fit_t(utils::TreeFitter fitter);

            static utils::TreeFitter Make(const ParticleTypeDatabase::Type& subtree);

            static TParticlePtr FindBest(const utils::TreeFitter::tree_t& fitted,
                                         const Particles_t& particles
                                         );

            utils::TreeFitter treefitter;
            utils::TreeFitter::tree_t fitted_Pi0;
            utils::TreeFitter::tree_t fitted_g1_Pi0;
            utils::TreeFitter::tree_t fitted_g2_Pi0;

            utils::TreeFitter::tree_t fitted_Omega;
            utils::TreeFitter::tree_t fitted_g_Omega;

            utils::TreeFitter::tree_t fitted_EtaPrime;
            utils::TreeFitter::tree_t fitted_g_EtaPrime;

        };

        struct Pi0_t : Fit_t {

            Pi0_t();

            struct BaseTree_t : Fit_t::BaseTree_t {
                ADD_BRANCH_T(std::vector<double>, IM_Pi0g, 2)
                void Reset();
            };

            BaseTree_t t;

            void Process(const Particles_t& particles, const TParticleTree_t& ptree_sigref);
        };

        struct OmegaPi0_t : Fit_t {

            OmegaPi0_t();

            struct BaseTree_t : Fit_t::BaseTree_t {
                ADD_BRANCH_T(double, IM_Pi0g_fitted)
                ADD_BRANCH_T(double, IM_Pi0g_best)

                ADD_BRANCH_T(double, Bachelor_E)

                void Reset();
            };

            BaseTree_t t;

            void Process(const Particles_t& particles, const TParticleTree_t& ptree_sigref);

        };

        Sig_t();

        Pi0_t Pi0;
        OmegaPi0_t OmegaPi0;

        struct SharedTree_t : EtapOmegaG::SharedTree_t {
            ADD_BRANCH_T(std::vector<double>, ggg, 4)
            ADD_BRANCH_T(std::vector<double>, gg_gg1, 3)
            ADD_BRANCH_T(std::vector<double>, gg_gg2, 3)

            ADD_BRANCH_T(double, AntiPi0FitProb)
            ADD_BRANCH_T(int,    AntiPi0FitIterations)

            ADD_BRANCH_T(double,   AntiPi0BeamEPull)
            ADD_BRANCH_T(double,   AntiPi0ProtonEPull)
            ADD_BRANCH_T(double,   AntiPi0ProtonThetaPull)
            ADD_BRANCH_T(double,   AntiPi0ProtonPhiPull)
            ADD_BRANCH_T(std::vector<double>,  AntiPi0PhotonEPulls)
            ADD_BRANCH_T(std::vector<double>,  AntiPi0PhotonThetaPulls)
            ADD_BRANCH_T(std::vector<double>,  AntiPi0PhotonPhiPulls)

            ADD_BRANCH_T(double, AntiEtaFitProb)
            ADD_BRANCH_T(int,    AntiEtaFitIterations)

            ADD_BRANCH_T(double,   AntiEtaBeamEPull)
            ADD_BRANCH_T(double,   AntiEtaProtonEPull)
            ADD_BRANCH_T(double,   AntiEtaProtonThetaPull)
            ADD_BRANCH_T(double,   AntiEtaProtonPhiPull)
            ADD_BRANCH_T(std::vector<double>,  AntiEtaPhotonEPulls)
            ADD_BRANCH_T(std::vector<double>,  AntiEtaPhotonThetaPulls)
            ADD_BRANCH_T(std::vector<double>,  AntiEtaPhotonPhiPulls)

            void Reset();
        };

        void SetupTrees(HistogramFactory HistFac);
        void Fill();
        void ResetBranches();
        void Process(const Particles_t& particles, const TParticleTree_t& ptree_sigref);

        SharedTree_t t;

    private:
        utils::TreeFitter treefitter_Pi0Pi0;
        utils::TreeFitter treefitter_Pi0Eta;
        void DoAntiPi0Eta(const Particles_t& particles);
        void DoPhotonCombinatorics(const TParticleList& photons);

    };

    struct Ref_t {

        struct Tree_t : EtapOmegaG::SharedTree_t {
            ADD_BRANCH_T(double,   IM_2g)

            void Reset();
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

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned max_iterations);

    EtapOmegaG(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
