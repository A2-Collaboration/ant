#pragma once

#include "analysis/utils/fitter/TreeFitter.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "utils/A2GeoAcceptance.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/TriggerSimulation.h"

#include "analysis/physics/scratch/wolfes/tools/tools.h"

#include "base/WrapTTree.h"

#include "TLorentzVector.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

struct singlePi0 :  Physics {

    //===================== Settings   ========================================================


    struct settings_t
    {
        const std::string Tree_Name = "tree";

        const unsigned nPhotons = 2;

        const interval<size_t>  Cut_NCands     = {3,3};
        const IntervalD         Cut_ProtonCopl = {-25,25};
        const IntervalD         Cut_MM         = ParticleTypeDatabase::Proton.GetWindow(350).Round();
        const IntervalD         Cut_MMAngle    = {0,25};
        const IntervalD         Cut_EMB_prob   = {0.005,1};

        const IntervalD              Range_Prompt  =   { -5,  5};
        const std::vector<IntervalD> Ranges_Random = { {-55,-10},
                                                       { 10, 55}  };

        const double fitter_ZVertex = 3;

        const unsigned Index_Data    = 0;
        const unsigned Index_Signal  = 1;
        const unsigned Index_MainBkg = 2;
        const unsigned Index_Offset  = 10;
        const unsigned Index_Unknown = 9;

    };

    const settings_t phSettings;
    const std::shared_ptr<TaggerDetector_t> tagger;


    ant::analysis::utils::A2SimpleGeometry geometry;

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
    static const std::vector<named_channel_t> otherBackgrounds;

    //===================== Histograms ========================================================

    TH1D* hist_steps             = nullptr;
    TH1D* hist_channels          = nullptr;
    TH1D* hist_channels_end      = nullptr;
    TH2D* hist_neutrals_channels = nullptr;

    //===================== KinFitting ========================================================


    std::shared_ptr<utils::UncertaintyModel> uncertModel = std::make_shared<utils::UncertaintyModels::FitterSergey>();

    utils::KinFitter kinFitterEMB;

    utils::TreeFitter fitterSig;
    std::vector<utils::TreeFitter::tree_t> pionsFitterSig;

    //========================  ProptRa. ============================================================

    utils::TriggerSimulation triggersimu;
    ant::analysis::PromptRandom::Switch promptrandom;

    //========================  Storage  ============================================================
    struct protonSelection_t
    {
        TParticlePtr   Proton;
        TParticleList  Photons;
        LorentzVec     PhotonSum;
        LorentzVec     Proton_MM;
        LorentzVec     PhotonBeam;
        double         Copl_pg;
        double         Angle_pMM;
        double         Tagg_E;

        template<typename wtf_ITER>
        protonSelection_t(const wtf_ITER& selectedProton, const TCandidateList& candidates,
                          const LorentzVec& photonBeam,   double taggE):
            PhotonSum({0,0,0},0),
            PhotonBeam(photonBeam),
            Tagg_E(taggE)
        {
            Proton = std::make_shared<TParticle>(ParticleTypeDatabase::Proton, selectedProton);
            for ( auto i_photon : candidates.get_iter())
                if (!(i_photon == selectedProton))
                {
                    Photons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));
                    PhotonSum += *Photons.back();
                }
            Proton_MM =   photonBeam
                        + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass())
                        - PhotonSum;
            Copl_pg   =   std_ext::radian_to_degree(vec2::Phi_mpi_pi(Proton->Phi()-PhotonSum.Phi() - M_PI ));
            Angle_pMM =   std_ext::radian_to_degree(Proton_MM.Angle(Proton->p));
        }
    };

    struct fitRatings_t
    {
        double Prob;
        double Chi2;
        int    Niter;
        bool   FitOk;
        std::vector<TLorentzVector> Intermediates;
        fitRatings_t(double prob,double chi2,int niter, bool fitOk,
                     const std::vector<TLorentzVector> intermediates):
            Prob(prob),Chi2(chi2),Niter(niter), FitOk(fitOk), Intermediates(intermediates){}
    };

    struct PionProdTree : WrapTTree
    {
        struct particleStorage_t
        {
            TParticlePtr   Proton;
            TParticleList  Photons;
            LorentzVec     PhotonSum;
            particleStorage_t():
                Proton(std::make_shared<TParticle>(ParticleTypeDatabase::Proton,LorentzVec({0,0,0},0))),
                Photons(TParticleList()),
                PhotonSum({0,0,0},0){}
            particleStorage_t(const TParticlePtr& proton,
                              const TParticleList& photons, const LorentzVec& photonSum):
                Proton(proton),
                Photons(photons),PhotonSum(photonSum){}
        };

        // type: 0   data
        //       1   signal (pi0)
        //       2   mainBkg(eta->gg)
        //       10+ otherBkg
        ADD_BRANCH_T(unsigned, MCTrue)

        ADD_BRANCH_T(double,   Tagg_W)
        ADD_BRANCH_T(unsigned, Tagg_Ch)
        ADD_BRANCH_T(double,   Tagg_E)
        ADD_BRANCH_T(double,   Tagg_Eff)
        ADD_BRANCH_T(double,   Tagg_EffErr)

        // sclowcontrol
        ADD_BRANCH_T(std::vector<double>,   TaggRates)
        ADD_BRANCH_T(double,                ExpLivetime)


        ADD_BRANCH_T(unsigned,   Neutrals)
//        ADD_BRANCH_T(double,   ChargedClusterE)
        ADD_BRANCH_T(double,   ProtonVetoE)
        ADD_BRANCH_T(double,   PionVetoE)

        ADD_BRANCH_T(double, CBAvgTime)
        ADD_BRANCH_T(double, CBESum)

        // best emb combination raw
        ADD_BRANCH_T(TLorentzVector,              proton)
        ADD_BRANCH_T(double,                      protonTime)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons)
        ADD_BRANCH_T(std::vector<double>,         photonTimes)
        ADD_BRANCH_T(TLorentzVector,              photonSum)
        ADD_BRANCH_T(double,                      IM2g)

        ADD_BRANCH_T(TLorentzVector,              proton_MM)
        ADD_BRANCH_T(double,                      pMM_angle)
        ADD_BRANCH_T(double,                      pg_copl)
        void SetRaw(const tools::protonSelection_t& selection);

        // best emb comb. emb-fitted
        ADD_BRANCH_T(TLorentzVector,              EMB_proton)
        ADD_BRANCH_T(std::vector<TLorentzVector>, EMB_photons)
        ADD_BRANCH_T(TLorentzVector,              EMB_photonSum)
        ADD_BRANCH_T(double,                      EMB_IM2g)
        ADD_BRANCH_T(double,                      EMB_Ebeam)

        ADD_BRANCH_T(double,                      EMB_prob)
        ADD_BRANCH_T(double,                      EMB_chi2)
        ADD_BRANCH_T(int,                         EMB_iterations)
        void SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result);

        //best tree-fit combination raw
        ADD_BRANCH_T(double,                        SIG_prob)
        ADD_BRANCH_T(double,                        SIG_chi2)
        ADD_BRANCH_T(int,                           SIG_iterations)
        ADD_BRANCH_T(std::vector<TLorentzVector>,   SIG_pions)
        void SetSIG(const singlePi0::fitRatings_t&  fitRating);

        static constexpr const char* treeName()       {return "tree";}
        static constexpr const char* treeAccessName() {return "singlePi0/tree";}
    };
    PionProdTree tree;

    //========================  MAIN     ============================================================

    singlePi0(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override {}
    virtual void ShowResult() override;

    //========================  TOOLS    ============================================================

    void FillStep(const std::string& step) {hist_steps->Fill(step.c_str(),1);}

};


}}}
