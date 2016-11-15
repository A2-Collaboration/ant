#pragma once

#include "analysis/utils/Fitter.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "utils/A2GeoAcceptance.h"

#include "base/WrapTTree.h"

#include "TLorentzVector.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

struct triplePi0 :  Physics {

    //===================== Settings   ========================================================


    struct settings_t
    {
        const std::string Tree_Name = "tree";

        const unsigned nPhotons = 6;

        const bool Opt_AllChannels;

        const double    Cut_CBESum     = 550;
        const unsigned  Cut_NCands     = 7;
        const IntervalD Cut_ProtonCopl = {-25,25};
        const IntervalD Cut_MM         = {850,1026};
        const double    Cut_MMAngle    = 20;
        const IntervalD Cut_EMB_Chi2    = {0.,40.};

        const IntervalD              Range_Prompt  =   { -5,  5};
        const std::vector<IntervalD> Ranges_Random = { {-55,-10},
                                                       { 10, 55}  };

        const double fitter_ZVertex = 3;

        const unsigned Index_Data    = 0;
        const unsigned Index_Signal  = 1;
        const unsigned Index_MainBkg = 2;
        const unsigned Index_Offset  = 10;
        const unsigned Index_Unknown = 9;


        settings_t(bool allChannels):
            Opt_AllChannels(allChannels){}
    };

    const settings_t phSettings;

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

    //===================== KinFitting ========================================================


    std::shared_ptr<utils::UncertaintyModel> uncertModel = std::make_shared<utils::UncertaintyModels::FitterSergey>();

    utils::KinFitter kinFitterEMB;

    utils::TreeFitter fitterSig;
    std::vector<utils::TreeFitter::tree_t> pionsFitterSig;

    utils::TreeFitter fitterBkg;
    std::vector<utils::TreeFitter::tree_t> pionsFitterBkg;


    //========================  ProptRa. ============================================================

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
        std::vector<TLorentzVector> Intermediates;
        fitRatings_t(double prob,double chi2,int niter,
                     const std::vector<TLorentzVector> intermediates):
            Prob(prob),Chi2(chi2),Niter(niter),Intermediates(intermediates){}
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
        //       1   signal (3pi0)
        //       2   mainBkg(eta->3pi0)
        //       10+ otherBkg
        ADD_BRANCH_T(unsigned, MCTrue)

        ADD_BRANCH_T(double,   Tagg_W)
        ADD_BRANCH_T(unsigned, Tagg_Ch)
        ADD_BRANCH_T(double,   Tagg_E)

        ADD_BRANCH_T(double, CBAvgTime)
        ADD_BRANCH_T(double, CBESum)

        // best emb combination raw
        ADD_BRANCH_T(TLorentzVector,              proton)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons)
        ADD_BRANCH_T(TLorentzVector,              photonSum)
        ADD_BRANCH_T(double,                      IM6g)

        ADD_BRANCH_T(TLorentzVector,              proton_MM)
        ADD_BRANCH_T(double,                      pMM_angle)
        ADD_BRANCH_T(double,                      pg_copl)
        void SetRaw(const triplePi0::protonSelection_t& selection);

        // best emb comb. emb-fitted
        ADD_BRANCH_T(TLorentzVector,              EMB_proton)
        ADD_BRANCH_T(std::vector<TLorentzVector>, EMB_photons)
        ADD_BRANCH_T(TLorentzVector,              EMB_photonSum)
        ADD_BRANCH_T(double,                      EMB_IM6g)
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
        void SetSIG(const triplePi0::fitRatings_t&  fitRating);

        ADD_BRANCH_T(double,                       BKG_prob)
        ADD_BRANCH_T(double,                       BKG_chi2)
        ADD_BRANCH_T(int,                          BKG_iterations)
        ADD_BRANCH_T(std::vector<TLorentzVector>,  BKG_pions)
        void SetBKG(const triplePi0::fitRatings_t& fitRating);
    };
    PionProdTree tree;

    //========================  MAIN     ============================================================

    triplePi0(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override {}
    virtual void ShowResult() override;

    //========================  TOOLS    ============================================================

    template<typename wtf_ITER>
    PionProdTree::particleStorage_t makeProtonSelection(const wtf_ITER& selectedProton, const TCandidateList& candidates)
    {
        const auto proton = std::make_shared<TParticle>(ParticleTypeDatabase::Proton, selectedProton);
        TParticleList photons;
        LorentzVec photonSum({0,0,0},0);
        for ( auto i_photon : candidates.get_iter())
            if (!(i_photon == selectedProton))
            {
                photons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));
                photonSum += *photons.back();
            }
        return PionProdTree::particleStorage_t(proton,photons,photonSum);
    }



    static std::vector<TLorentzVector> MakeTLorenz(const TParticleList& particles)
    {
        std::vector<TLorentzVector> lg(particles.size());
        std::transform(particles.begin(),particles.end(),lg.begin(),
                       [](const TParticlePtr& ph){return TLorentzVector(*ph);});
        return lg;
    }

    void FillStep(const std::string& step) {hist_steps->Fill(step.c_str(),1);}

};


}}}
