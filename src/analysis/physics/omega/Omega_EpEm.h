#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "base/WrapTTree.h"
#include "utils/TriggerSimulation.h"

namespace ant {
namespace analysis {
namespace physics {

class Omega_EpEm : public Physics
{
public:
    struct tree_t : WrapTTree {
        ADD_BRANCH_T(bool,      IsMC)
        ADD_BRANCH_T(double,    TaggW)
        ADD_BRANCH_T(unsigned,  nClusters)
    };


public:
    Omega_EpEm(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

    virtual void ShowResult() override;

protected:

    TH1D* h_nClusters;
    TH1D* h_nClusters_pr;
    TH1D* h_nPhotons;
    TH1D* h_IM_2g;
    tree_t t;
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    int     nCB = 0;
    int     nTAPS = 0;
    double  CBAvgTime  = 0.0;
    double  CBSumVetoE  = 0.0;

};


}}} // end of ant::analysis::physics
