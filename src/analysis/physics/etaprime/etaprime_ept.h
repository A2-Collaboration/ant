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

struct EtapEPT : Physics {

    TH1D* h_CommonCuts;

    // TreeCommon contains things
    // shared among sig/ref analyses

    struct TreeCommon : WrapTTree {
        ADD_BRANCH_T(unsigned, MCTrue)
        ADD_BRANCH_T(double,   TrueZVertex)
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

        ADD_BRANCH_T(std::vector<double>,  PhotonThetas)

        ADD_BRANCH_T(double,   KinFitProb)
        ADD_BRANCH_T(int,      KinFitIterations)
        ADD_BRANCH_T(double,   KinFitZVertex)

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

    utils::KinFitter kinfitter_ref;

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
                              Particles_t& particles,
                              SharedTree_t& t,
                              TH1D* h_CommonCuts
                              );

    static bool doKinfit(const TTaggerHit& taggerhit,
                         utils::KinFitter& kinfitter,
                         Particles_t& particles,
                         SharedTree_t& t,
                         TH1D* h_CommonCuts);

     struct Ref_t {

        struct Tree_t : EtapEPT::SharedTree_t {
            ADD_BRANCH_T(double,   IM_2g)

            void Reset();
        };
        Tree_t t;

        void ResetBranches();
        void Process(const Particles_t& particles);
    };


    Ref_t Ref;

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned max_iterations);

    EtapEPT(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
