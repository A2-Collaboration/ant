#pragma once

#include "analysis/utils/fitter/TreeFitter.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "utils/A2GeoAcceptance.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/TriggerSimulation.h"

#include "utils/ProtonPhotonCombs.h"

#include "analysis/physics/scratch/wolfes/tools/tools.h"

#include "base/WrapTTree.h"

#include "TLorentzVector.h"
#include "tree/TSimpleParticle.h"

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

        const interval<size_t>  Cut_NCands         = {3,10};
        const IntervalD         Cut_ProtonCopl     = {-25,25};
        const IntervalD         Cut_MM             = ParticleTypeDatabase::Proton.GetWindow(350).Round();
        const IntervalD         Cut_IM             = ParticleTypeDatabase::Pi0.GetWindow(ParticleTypeDatabase::Pi0.Mass()).Round();
        const IntervalD         Cut_MMAngle        = {0,25};
        const IntervalD         Cut_EMB_prob       = {0.005,1};

        const double            Cut_MaxDiscardedEk = 100.;


        const IntervalD              Range_Prompt  =   { -5,  5};
        const std::vector<IntervalD> Ranges_Random = { {-55,-10},
                                                       { 10, 55}  };
        const double fitter_ZVertex = 3;

        const unsigned Index_Data       = 0;
        const unsigned Index_Signal     = 1;
        const unsigned Index_MainBkg    = 2;
        const unsigned Index_Offset     = 11;
        const unsigned Index_UnTagged   = 10;
        const unsigned Index_brokenTree = 9;

    };

    bool FinalCuts() const { return tree.Neutrals() < 2;}

    const settings_t phSettings;
    const bool flag_mc;
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
    TH1D* hist_tagger_hits       = nullptr;

    TH1D* hist_ncands            = nullptr;

    TH1D* seenMC              = nullptr;
    TH1D* taggerScalars       = nullptr;


    TH2D* hist_seen       = nullptr;
    TH2D* hist_rec        = nullptr;
    TH2D* hist_efficiency = nullptr;

    //===================== KinFitting ========================================================


    std::shared_ptr<utils::UncertaintyModel> uncertModelData = std::make_shared<utils::UncertaintyModels::FitterSergey>();
    std::shared_ptr<utils::UncertaintyModel> uncertModelMC = std::make_shared<utils::UncertaintyModels::FitterSergey>();

    utils::KinFitter fitterEMB;


    //========================  ProptRa. ============================================================
    utils::TriggerSimulation triggersimu;
    ant::analysis::PromptRandom::Switch promptrandom;


    struct effTree_t : WrapTTree
    {
        ADD_BRANCH_T(double, Theta)
        ADD_BRANCH_T(double, Phi)
        ADD_BRANCH_T(double, CosThetaPi0)
        ADD_BRANCH_T(double, Egamma)
        ADD_BRANCH_T(int,    TaggerBin)

        virtual std::string treeName() const=0;
        virtual std::string treeAccessName() const {return "singlePi0/" + treeName();}

        virtual ~effTree_t(){}

    };

    struct SeenTree : effTree_t
    {
        virtual std::string treeName() const override {return "seen";}
    };
    SeenTree seenSignal;

    struct RecTree : effTree_t
    {
        virtual std::string treeName() const override {return "rec";}
    };
    RecTree recSignal;



    struct PionProdTree : WrapTTree
    {
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
        ADD_BRANCH_T(double,                ExpLivetime)


        ADD_BRANCH_T(unsigned,   Neutrals)
        ADD_BRANCH_T(double,     PionVetoE)

        ADD_BRANCH_T(double,     CBAvgTime)
        ADD_BRANCH_T(double,     CBESum)

        ADD_BRANCH_T(int,        NCands)

        // best emb combination raw
        ADD_BRANCH_T(TSimpleParticle,              proton)
        ADD_BRANCH_T(std::vector<TSimpleParticle>, photons)

        ADD_BRANCH_T(TLorentzVector,               photonSum)
        ADD_BRANCH_T(double,                       IM2g)
        ADD_BRANCH_T(double,                       cosThetaPi0COMS)

        ADD_BRANCH_T(double,                      IMproton_MM)
        ADD_BRANCH_T(double,                      DiscardedEk)
        ADD_BRANCH_T(TLorentzVector,              proton_MM)
        ADD_BRANCH_T(double,                      pMM_angle)
        ADD_BRANCH_T(double,                      pg_copl)
        void SetRaw(const utils::ProtonPhotonCombs::comb_t& selection);

        // best emb comb. emb-fitted
        ADD_BRANCH_T(TSimpleParticle,              EMB_proton)
        ADD_BRANCH_T(std::vector<TSimpleParticle>, EMB_photons)
        ADD_BRANCH_T(TLorentzVector,               EMB_photonSum)
        ADD_BRANCH_T(double,                       EMB_IM2g)
        ADD_BRANCH_T(double,                       EMB_Ebeam)
        ADD_BRANCH_T(double,                       EMB_cosThetaPi0COMS)

        ADD_BRANCH_T(double,                       EMB_prob)
        ADD_BRANCH_T(double,                       EMB_chi2)
        ADD_BRANCH_T(int,                          EMB_iterations)
        void SetEMB(const utils::KinFitter& kF, const APLCON::Result_t& result);

        static constexpr const char* treeName()       {return "tree";}
        static constexpr const char* treeAccessName() {return "singlePi0/tree";}
    };
    PionProdTree tree;

    //========================  MAIN     ============================================================

    singlePi0(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

    //========================  TOOLS    ============================================================

    void FillStep(const std::string& step) {hist_steps->Fill(step.c_str(),1);}

};


}}}
