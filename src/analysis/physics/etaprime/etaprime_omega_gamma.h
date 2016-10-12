#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/utils/Fitter.h"
#include "analysis/utils/MCSmear.h"
#include "analysis/utils/A2GeoAcceptance.h"
#include "analysis/plot/PromptRandomHist.h"

#include "base/ParticleTypeTree.h"
#include "base/WrapTTree.h"

#include "TLorentzVector.h"

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

    TH1D* h_Cuts;

    TH1D* h_LostPhotons_sig;
    TH1D* h_LostPhotons_ref;
    TH1D* h_MissedBkg;

    struct TreeCommon : virtual WrapTTree {
        ADD_BRANCH_T(unsigned, MCTrue)
        ADD_BRANCH_T(double,   TrueZVertex)
        ADD_BRANCH_T(double,   CBSumE)
        ADD_BRANCH_T(double,   CBAvgTime)
        ADD_BRANCH_T(double,   PIDSumE)

        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)
    };

    PromptRandom::Switch promptrandom;

    struct fitparams_t {
        const utils::UncertaintyModelPtr Fit_uncertainty_model;
        const bool Fit_Z_vertex;
        const double Z_vertex_sigma;
        fitparams_t(utils::UncertaintyModelPtr fit_uncertainty_model,
                    bool fit_Z_vertex, double Z_vertex_sigma) :
            Fit_uncertainty_model(fit_uncertainty_model),
            Fit_Z_vertex(fit_Z_vertex),
            Z_vertex_sigma(Z_vertex_sigma)
        {}
    };

    const bool FlagWolfgang;
    const fitparams_t fitparams;

    std::unique_ptr<utils::MCSmear> mc_smear;

    utils::A2SimpleGeometry geometry;

    // TreeCommon contains things
    // shared among sig/ref analyses
    TreeCommon t;

    struct particle_t {
        particle_t(const TParticlePtr& proton) : Proton(proton) {}
        TParticleList Photons;
        TParticlePtr  Proton;

        double     MissingMass = std_ext::NaN;
        LorentzVec PhotonSum{};
        double     DiscardedEk = 0;
    };

    struct params_t {
        bool Filter(
                unsigned n, TH1D* h_Cuts,
                double maxDiscardedEk = std_ext::inf,
                interval<double> missingMassCut = {-std_ext::inf, std_ext::inf},
                interval<double> photonSumCut = {-std_ext::inf, std_ext::inf}
                );

        std::list<particle_t> Particles{};

        TTaggerHit TaggerHit{};

        TParticleTree_t ParticleTree = nullptr;
        bool IsSignalTree = false;
    };

    struct ProtonPhotonTree_t : virtual WrapTTree {
        ADD_BRANCH_T(double,   DiscardedEk)

        ADD_BRANCH_T(double,   PhotonsEk)
        ADD_BRANCH_T(unsigned, nPhotonsCB)
        ADD_BRANCH_T(unsigned, nPhotonsTAPS)
        ADD_BRANCH_T(double,   CBSumVetoE)
        ADD_BRANCH_T(double,   PhotonSum)
        ADD_BRANCH_T(std::vector<double>,  PhotonThetas)

        ADD_BRANCH_T(double,   ProtonTime)
        ADD_BRANCH_T(double,   ProtonE)
        ADD_BRANCH_T(double,   ProtonTheta)
        ADD_BRANCH_T(double,   ProtonVetoE)
        ADD_BRANCH_T(double,   ProtonShortE)
        ADD_BRANCH_T(double,   ProtonTrueAngle)
        ADD_BRANCH_T(double,   ProtonCopl)
        ADD_BRANCH_T(double,   MissingMass)

        ADD_BRANCH_T(double,   FittedProtonE)

        void Fill(const params_t& params, const particle_t& p, double fitted_proton_E);

    };

    struct Sig_t {

        // the subtree to be fitted is either
        // pi0->2g or omega->pi0g->3g
        // fitting the whole decay tree would overconstrain the
        // photons

        struct SharedTree_t : virtual WrapTTree {
            ADD_BRANCH_T(double, KinFitProb)
            ADD_BRANCH_T(int,    KinFitIterations)
            ADD_BRANCH_T(double, KinFitZVertex)

            ADD_BRANCH_T(double, AntiPi0FitProb)
            ADD_BRANCH_T(int,    AntiPi0FitIterations)
            ADD_BRANCH_T(double, AntiPi0FitZVertex)

            ADD_BRANCH_T(double, AntiEtaFitProb)
            ADD_BRANCH_T(int,    AntiEtaFitIterations)
            ADD_BRANCH_T(double, AntiEtaFitZVertex)

        };

        struct Fit_t {

            struct BaseTree_t : virtual ProtonPhotonTree_t {

                ADD_BRANCH_T(std::vector<double>, ggg, 4)
                ADD_BRANCH_T(std::vector<double>, gg_gg1, 3)
                ADD_BRANCH_T(std::vector<double>, gg_gg2, 3)

                ADD_BRANCH_T(double, TreeFitProb)
                ADD_BRANCH_T(int,    TreeFitIterations)
                ADD_BRANCH_T(double, TreeFitZVertex)

                ADD_BRANCH_T(double, IM_Pi0)
                ADD_BRANCH_T(double, IM_Pi0gg)
                ADD_BRANCH_T(double, IM_gg)

                ADD_BRANCH_T(unsigned, MCTrueMatch)

                // information about the two photons NOT belonging to the Pi0
                ADD_BRANCH_T(std::vector<double>,   gNonPi0_Theta, 2)
                ADD_BRANCH_T(std::vector<double>,   gNonPi0_CaloE, 2)

            };

            Fit_t(utils::TreeFitter fitter);

            void FillPhotonCombs(BaseTree_t& t, const TParticleList& photons);

            static utils::TreeFitter Make(const ParticleTypeDatabase::Type& subtree, fitparams_t fitparams);

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

            const bool FlagWolfgang;

            Pi0_t(fitparams_t fitparams, bool flagWolfgang);

            struct Tree_t : virtual Fit_t::BaseTree_t {
                ADD_BRANCH_T(std::vector<double>, IM_Pi0g, 2)
                ADD_BRANCH_T(std::vector<double>, Bachelor_E, 2)
            };

            Tree_t t;

            struct WolfgangsTree_t : TreeCommon, SharedTree_t, Pi0_t::Tree_t {
                ADD_BRANCH_T(TLorentzVector, Proton)
                ADD_BRANCH_T(TLorentzVector, Photon1)
                ADD_BRANCH_T(TLorentzVector, Photon2)
                ADD_BRANCH_T(TLorentzVector, Photon3)
                ADD_BRANCH_T(TLorentzVector, Photon4)
            };

            WolfgangsTree_t t_w;

            void Process(const params_t& params);
        };

        struct OmegaPi0_t : Fit_t {

            OmegaPi0_t(fitparams_t fitparams);

            struct Tree_t : Fit_t::BaseTree_t {
                ADD_BRANCH_T(double, IM_Pi0g)
                ADD_BRANCH_T(double, Bachelor_E)

            };

            Tree_t t;

            void Process(const params_t& p);

        };

        const bool FlagWolfgang;

        Sig_t(const HistogramFactory& HistFac, fitparams_t fitparams, bool flagWolfgang);

        TH1D* h_Cuts;

        TTree* treeCommon;
        SharedTree_t t;

        Pi0_t Pi0;
        OmegaPi0_t OmegaPi0;

        utils::KinFitter  kinfitter;
        utils::TreeFitter treefitter_Pi0Pi0;
        utils::TreeFitter treefitter_Pi0Eta;

        void Process(params_t params);
        void DoAntiPi0Eta(const params_t& params);

    };

    struct Ref_t {

        struct Tree_t : EtapOmegaG::ProtonPhotonTree_t {
            ADD_BRANCH_T(double, KinFitProb)
            ADD_BRANCH_T(int,    KinFitIterations)
            ADD_BRANCH_T(double, KinFitZVertex)

            ADD_BRANCH_T(double,   IM_2g)
        };

        Ref_t(const HistogramFactory& HistFac, fitparams_t fitparams, bool flagWolfgang);

        TH1D* h_Cuts;

        TTree* treeCommon;
        Tree_t t;

        utils::KinFitter kinfitter;

        void Process(params_t params);
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
