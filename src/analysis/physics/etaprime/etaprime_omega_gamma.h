#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/utils/Fitter.h"
#include "analysis/utils/MCSmear.h"
#include "analysis/utils/A2GeoAcceptance.h"
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

    TH1D* h_LostPhotons_sig;
    TH1D* h_LostPhotons_ref;
    TH1D* h_MissedBkg;

    // TreeCommon contains things
    // shared among sig/ref analyses

    struct TreeCommon : WrapTTree {
        ADD_BRANCH_T(unsigned, MCTrue)
        ADD_BRANCH_T(double,   TrueZVertex)
        ADD_BRANCH_T(double,   CBSumE)
        ADD_BRANCH_T(double,   CBAvgTime)
        ADD_BRANCH_T(double,   PIDSumE)

        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)

        ADD_BRANCH_T(unsigned, TaggNPrompt)
        ADD_BRANCH_T(unsigned, TaggNRandom)
    };

    struct SharedTree_t : WrapTTree {
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

        ADD_BRANCH_T(std::vector<double>,  PhotonThetas)

        ADD_BRANCH_T(double,   KinFitProb)
        ADD_BRANCH_T(int,      KinFitIterations)
        ADD_BRANCH_T(double,   KinFitZVertex)

        ADD_BRANCH_T(double, KinFitBeamEPull)
        ADD_BRANCH_T(std::vector<double>, KinFitProtonPulls)
        ADD_BRANCH_T(std::vector<std::vector<double>>, KinFitPhotonsPulls)
    };

    TreeCommon t;

    PromptRandom::Switch promptrandom;

    struct params_t {
        const utils::UncertaintyModelPtr Fit_uncertainty_model;
        const bool Fit_Z_vertex;
        const double Z_vertex_sigma;
        params_t(utils::UncertaintyModelPtr fit_uncertainty_model,
                 bool fit_Z_vertex, double Z_vertex_sigma) :
            Fit_uncertainty_model(fit_uncertainty_model),
            Fit_Z_vertex(fit_Z_vertex),
            Z_vertex_sigma(Z_vertex_sigma)
        {}
    };

    const params_t params;

    utils::KinFitter kinfitter_sig;
    utils::KinFitter kinfitter_ref;

    std::unique_ptr<utils::MCSmear>             mc_smear;
    std::unique_ptr<utils::MCFakeReconstructed> mc_fake;

    utils::A2SimpleGeometry geometry;

    struct Particles_t {
        double         PhotonEnergy;
        TParticlePtr   Proton;
        TParticleList  Photons;
        TParticleList  FittedPhotons;
        LorentzVec     PhotonSum;
        LorentzVec     FittedPhotonSum;
        double         DiscardedEk = 0;
    };

    static bool doKinfit(const TTaggerHit& taggerhit,
                         TParticlePtr true_proton,
                         utils::KinFitter& kinfitter,
                         Particles_t& particles,
                         SharedTree_t& t,
                         TH1D* h_CommonCuts);

    struct Sig_t {

        // the subtree to be fitted is either pi0->2g
        // or omega->pi0g->3g
        // fitting the whole decay tree would overconstrain the
        // photons

        struct Fit_t {

            struct BaseTree_t : WrapTTree {

                ADD_BRANCH_T(double, TreeFitProb)
                ADD_BRANCH_T(int,    TreeFitIterations)
                ADD_BRANCH_T(double, TreeFitZVertex)

                ADD_BRANCH_T(double, TreeFitBeamEPull)
                ADD_BRANCH_T(std::vector<double>, TreeFitProtonPulls)
                ADD_BRANCH_T(std::vector<std::vector<double>>, TreeFitPhotonsPulls)

                ADD_BRANCH_T(double, IM_Pi0_fitted)
                ADD_BRANCH_T(double, IM_Pi0_best)

                ADD_BRANCH_T(double, IM_Pi0gg_fitted)
                ADD_BRANCH_T(double, IM_Pi0gg_best)

                ADD_BRANCH_T(double, IM_gg_fitted)
                ADD_BRANCH_T(double, IM_gg_best)

                ADD_BRANCH_T(unsigned, MCTrueMatch)

                void Reset();
            };

        protected:

            Fit_t(utils::TreeFitter fitter);

            static utils::TreeFitter Make(const ParticleTypeDatabase::Type& subtree, params_t params);

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

            Pi0_t(params_t params);

            struct BaseTree_t : Fit_t::BaseTree_t {
                ADD_BRANCH_T(std::vector<double>, IM_Pi0g_fitted, 2)
                ADD_BRANCH_T(std::vector<double>, IM_Pi0g_best, 2)

                ADD_BRANCH_T(std::vector<double>, Bachelor_E_fitted, 2)
                ADD_BRANCH_T(std::vector<double>, Bachelor_E_best, 2)

                void Reset();
            };

            BaseTree_t t;

            void Process(const Particles_t& particles, const TParticleTree_t& ptree_sig);
        };

        struct OmegaPi0_t : Fit_t {

            OmegaPi0_t(params_t params);

            struct BaseTree_t : Fit_t::BaseTree_t {
                ADD_BRANCH_T(double, IM_Pi0g_fitted)
                ADD_BRANCH_T(double, IM_Pi0g_best)

                ADD_BRANCH_T(double, Bachelor_E_fitted)
                ADD_BRANCH_T(double, Bachelor_E_best)

                void Reset();
            };

            BaseTree_t t;

            void Process(const Particles_t& particles, const TParticleTree_t& ptree_sig);

        };

        Sig_t(params_t params);

        Pi0_t Pi0;
        OmegaPi0_t OmegaPi0;

        struct SharedTree_t : EtapOmegaG::SharedTree_t {
            ADD_BRANCH_T(std::vector<double>, ggg, 4)
            ADD_BRANCH_T(std::vector<double>, gg_gg1, 3)
            ADD_BRANCH_T(std::vector<double>, gg_gg2, 3)

            ADD_BRANCH_T(double, AntiPi0FitProb)
            ADD_BRANCH_T(int,    AntiPi0FitIterations)
            ADD_BRANCH_T(double, AntiPi0FitZVertex)

            ADD_BRANCH_T(double,   AntiPi0BeamEPull)
            ADD_BRANCH_T(std::vector<double>, AntiPi0ProtonPulls)
            ADD_BRANCH_T(std::vector<std::vector<double>>, AntiPi0PhotonsPulls)

            ADD_BRANCH_T(double, AntiEtaFitProb)
            ADD_BRANCH_T(int,    AntiEtaFitIterations)
            ADD_BRANCH_T(double, AntiEtaFitZVertex)


            ADD_BRANCH_T(double,   AntiEtaBeamEPull)
            ADD_BRANCH_T(std::vector<double>, AntiEtaProtonPulls)
            ADD_BRANCH_T(std::vector<std::vector<double>>, AntiEtaPhotonsPulls)

            void Reset();
        };

        void SetupTrees(HistogramFactory HistFac);
        void Fill();
        void ResetBranches();
        void Process(const Particles_t& particles, const TParticleTree_t& ptree_sig);

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
