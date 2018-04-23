#pragma once

#include <vector>

#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/TriggerSimulation.h"
#include "base/WrapTTree.h"
#include "tree/TSimpleParticle.h"

#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace physics {

class Omega_EpEm : public Physics {
public:
    struct tree_t : WrapTTree {
        ADD_BRANCH_T(bool,          IsMC)
        ADD_BRANCH_T(double,        TaggW)
        ADD_BRANCH_T(short,      nClusters)
        ADD_BRANCH_T(int,           nTAPSneutral)
        ADD_BRANCH_T(int,           nTAPScharged)
        ADD_BRANCH_T(int,           nCBneutral)
        ADD_BRANCH_T(int,           nCBcharged)


        ADD_BRANCH_T(std::vector<TSimpleParticle>,   p_tapsCharged)
        ADD_BRANCH_T(std::vector<TSimpleParticle>,   p_cbCharged)

        void fillAndReset() // to fill the tree and reset all values which add up
        {
            Tree->Fill();
            nTAPSneutral    = -1;
            nTAPScharged    = -1;
            nCBneutral      = -1;
            nCBcharged      = -1;
//            IsMC            = false;
//            nClusters       = 0;
//            TaggW           = 0;
//            p_tapsCharged().resize(0);
//            p_cbCharged().resize(0);
        }
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
    TH1D* h_nCandCBcharged = nullptr;
    TH1D* h_nCandTAPScharged = nullptr;
    TH1D* h_nClusters_pr = nullptr;
    TH1D* h_PIDenergy = nullptr;
    TH1D* h_TAPSVetoEnergy = nullptr;
    TH1D* h_IM = nullptr;
    TH1D* energy = nullptr;
    TH1D* theta = nullptr;
    TH1D* phi = nullptr;
    TH1D* detectors = nullptr;

    TH2D* h_nCandCharged = nullptr;
    TH2D* h_nCand = nullptr;
    TH2D* h_cbdEE = nullptr;
    tree_t t;
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    double  CBAvgTime  = 0.0;
    double  CBSumVetoE  = 0.0;


};


}}} // end of ant::analysis::physics
