#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "base/WrapTTree.h"

namespace ant {
namespace analysis {
namespace physics {

class Omega_EpEm : public Physics
{
public:
    struct tree_t : WrapTTree {
        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(unsigned, nClusters)
    };


public:
    Omega_EpEm(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

    virtual void ShowResult() override;

protected:

    TH1D* h_nClusters;
    TH1D* h_nClusters_pr;
    tree_t t;
    PromptRandom::Switch promptrandom;
    int     b_nCB = 0;
    int     b_nTAPS = 0;
    double  b_CBAvgTime  = 0.0;
    double  b_CBSumVetoE  = 0.0;

};


}}} // end of ant::analysis::physics
