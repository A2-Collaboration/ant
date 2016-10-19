#pragma once

#include "analysis/utils/Fitter.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "utils/A2GeoAcceptance.h"

#include "base/WrapTTree.h"


namespace ant {
namespace analysis {
namespace physics {

struct triplePi0 :  Physics {

    struct settings_t
    {
        const std::map<int,std::string> EventTypes= {{0,"signal"},        // 3 pi0 photoproduction
                                                     {1,"background"},    // eta -> pi0 pi0 pi0
                                                     {-1,"other"}};
        const double    CBESum        = 550;
        const unsigned  NCands        = 7;

        const IntervalD Cut_ProtonCoplCut = {-19,19};
        const double    Cut_MMAngleCut    = 15;
        const IntervalD Cut_EMBProbP   = {0.9,1};

        const IntervalD              Range_Prompt  =   { -5,  5};
        const std::vector<IntervalD> Ranges_Random = { {-55,-10},
                                                       { 10, 55}  };

        const double fitter_ZVertex = 3;

    };
    const settings_t phSettings;

    //geometry
    ant::analysis::utils::A2SimpleGeometry geometry;

    // =======================   constants =====================================================

    ParticleTypeTree signal_tree;
    ParticleTypeTree bkg_tree;



    //===================== KinFitting ========================================================

    std::shared_ptr<utils::UncertaintyModel> uncertModel = std::make_shared<utils::UncertaintyModels::FitterSergey>();

    utils::TreeFitter fitterSig;
    std::vector<utils::TreeFitter::tree_t> intermediatesTreeSig= std::vector<utils::TreeFitter::tree_t>(3);

    utils::TreeFitter fitterBkg;
    std::vector<utils::TreeFitter::tree_t> intermediatesTreeBkg= std::vector<utils::TreeFitter::tree_t>(3);

    utils::KinFitter kinFitterEMB;

    ant::analysis::PromptRandom::Switch promptrandom;



    //========================  Storage  ============================================================

    struct PionProdTree : WrapTTree
    {
        ADD_BRANCH_T(bool, isMC)

        ADD_BRANCH_T(double,   Tagg_W)
        ADD_BRANCH_T(unsigned, Tagg_Ch)
        ADD_BRANCH_T(double,   Tagg_E)

        ADD_BRANCH_T(double, CBAvgTime)

        ADD_BRANCH_T(std::vector<double>, ggIM)
    };




    triplePi0(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override {}
    virtual void ShowResult() override {}

};


}}}
