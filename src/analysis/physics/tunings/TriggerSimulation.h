#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/fitter/KinFitter.h"
#include "utils/TriggerSimulation.h"
#include "base/WrapTTree.h"

namespace ant {
namespace analysis {
namespace physics {

class TriggerSimulation : public Physics {
public:

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(bool,   IsMC)
        ADD_BRANCH_T(bool,   Triggered)
        ADD_BRANCH_T(double, CBEnergySum)

        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(unsigned, TaggCh)

        ADD_BRANCH_T(unsigned, nPhotons)
        ADD_BRANCH_T(double,   FitProb)
        ADD_BRANCH_T(double,   ZVertex)

        ADD_BRANCH_T(std::vector<double>, IM_Combs_fitted)
        ADD_BRANCH_T(std::vector<double>, IM_Combs_raw)

    };

protected:
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    TH1D* steps;

    TH1D* h_CBESum_raw;
    TH1D* h_CBESum_pr;
    TH1D* h_CBESum_fit;
    TH1D* h_CBTiming;
    TH2D* h_CBTiming_CaloE;


    TH1D* h_TaggT;
    TH1D* h_TaggT_corr;
    TH2D* h_TaggT_CBTiming;

    Tree_t t;

    utils::KinFitter fitter;

public:
    TriggerSimulation(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}}