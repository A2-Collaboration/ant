#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "base/WrapTTree.h"
#include "utils/TriggerSimulation.h"
#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace physics {

class Omega_EpEm : public Physics
{
public:
    struct tree_t : WrapTTree {
        ADD_BRANCH_T(bool,          IsMC)
        ADD_BRANCH_T(double,        TaggW)
        ADD_BRANCH_T(unsigned,      nClusters)


//        ADD_BRANCH_T(std::vector<TCandidate>,   p_tapsCharged)
//        ADD_BRANCH_T(std::vector<TCandidate>,   p_cbCharged)
    };


public:
    Omega_EpEm(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

    virtual void ShowResult() override;

protected:

    TH1D* h_nClusters = nullptr;
    TH1D* h_nCandidatesEvent = nullptr;
    TH1D* h_nCandCB = nullptr;
    TH1D* h_nCandTAPS = nullptr;
    TH1D* h_nClusters_pr = nullptr;
    TH1D* h_PIDenergy = nullptr;
    TH1D* h_TAPSVetoEnergy = nullptr;
    TH1D* h_IM = nullptr;
    TH1D* energy = nullptr;
    TH1D* theta = nullptr;
    TH1D* phi = nullptr;
    TH1D* detectors = nullptr;

    TH2D* h_cbdEE = nullptr;
    tree_t t;
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    double  CBAvgTime  = 0.0;
    double  CBSumVetoE  = 0.0;


};


}}} // end of ant::analysis::physics
