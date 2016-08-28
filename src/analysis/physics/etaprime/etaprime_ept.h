#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/Fitter.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/MCSmear.h"

#include "base/WrapTTree.h"

class TH1D;
class TH2D;
class TH3D;

namespace ant {

namespace expconfig {
namespace detector {
struct EPT;
}
}

namespace analysis {
namespace physics {

struct EtapEPT : Physics {

    std::shared_ptr<expconfig::detector::EPT> EPT;

    TH1D* h_Cuts;

    // TreeCommon contains things
    // shared among sig/ref analyses

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(double,   CBSumE)
        ADD_BRANCH_T(double,   CBAvgTime)
        ADD_BRANCH_T(double,   PIDSumE)

        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggW_wide)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(double,   TaggE_)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)
        ADD_BRANCH_T(unsigned, TaggCh_)


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

        ADD_BRANCH_T(double,   PhotonSum)
        ADD_BRANCH_T(double,   ProtonCopl)
        ADD_BRANCH_T(double,   MissingMass)
        ADD_BRANCH_T(double,   FittedProtonE)

        ADD_BRANCH_T(std::vector<double>,  PhotonThetas)

        ADD_BRANCH_T(double,   KinFitProb)
        ADD_BRANCH_T(int,      KinFitIterations)
        ADD_BRANCH_T(double,   KinFitZVertex)

        ADD_BRANCH_T(double,   KinFitBeamEPull)
        ADD_BRANCH_T(std::vector<double>, KinFitProtonPulls)
        ADD_BRANCH_T(std::vector<std::vector<double>>, KinFitPhotonsPulls)

        ADD_BRANCH_T(double,   IM_2g)
    };

    const bool taggChPerm;

    Tree_t t;

    PromptRandom::Switch promptrandom;
    PromptRandom::Switch promptrandom_wide;

    utils::KinFitter kinfitter;

    std::unique_ptr<utils::MCSmear> mc_smear;


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
                              Tree_t& t,
                              TH1D* h_Cuts
                              );

    static bool doKinfit(const TTaggerHit& taggerhit,
                         utils::KinFitter& kinfitter,
                         Particles_t& particles,
                         Tree_t& t,
                         TH1D* h_Cuts);

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned max_iterations);

    EtapEPT(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
