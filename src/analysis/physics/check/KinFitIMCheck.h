#pragma once

#include "analysis/physics/Physics.h"
#include "utils/TriggerSimulation.h"
#include "plot/PromptRandomHist.h"
#include "base/WrapTTree.h"
#include "utils/fitter/KinFitter.h"
#include <string>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class KinFitIMCheck : public Physics {
protected:
    utils::TriggerSimulation triggersimu;
    PromptRandom::Switch promptrandom;

    const utils::UncertaintyModelPtr fitmodel_data;
    const utils::UncertaintyModelPtr fitmodel_mc;
    utils::KinFitter fitter;

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(unsigned, nCB)
        ADD_BRANCH_T(unsigned, nTAPS)
        ADD_BRANCH_T(double, PIDSumE)
        ADD_BRANCH_T(double, CBVetoSumE)

        ADD_BRANCH_T(double, TaggW)
        ADD_BRANCH_T(unsigned, TaggCh)

        ADD_BRANCH_T(double, FitProb)
        ADD_BRANCH_T(double, ZVertex)

        ADD_BRANCH_T(double, IM_Raw)
        ADD_BRANCH_T(std::vector<double>, IM_Raw_1)
        ADD_BRANCH_T(std::vector<double>, IM_Raw_2)
        ADD_BRANCH_T(std::vector<double>, IM_Raw_3)
        ADD_BRANCH_T(std::vector<double>, IM_Raw_4)

        ADD_BRANCH_T(double, IM)
        ADD_BRANCH_T(std::vector<double>, IM_1)
        ADD_BRANCH_T(std::vector<double>, IM_2)
        ADD_BRANCH_T(std::vector<double>, IM_3)
        ADD_BRANCH_T(std::vector<double>, IM_4)
    };

    Tree_t t;

    TH1D* h_Steps;

public:
    KinFitIMCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t&) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics