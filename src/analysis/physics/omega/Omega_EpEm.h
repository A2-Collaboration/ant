#pragma once

#include <vector>

#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/ClusterTools.h"
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
        ADD_BRANCH_T(short,         nClusters)
        ADD_BRANCH_T(int,           nTAPSneutral)
        ADD_BRANCH_T(int,           nTAPScharged)
        ADD_BRANCH_T(int,           nCBneutral)
        ADD_BRANCH_T(int,           nCBcharged)
        ADD_BRANCH_T(std::vector<double>,   CB_effectiveradius)
        ADD_BRANCH_T(std::vector<double>,   CB_lateralmoment)
        ADD_BRANCH_T(int,           nCombsInIM)
        ADD_BRANCH_T(std::vector<TSimpleParticle>,   p_tapsCharged)
        ADD_BRANCH_T(std::vector<TSimpleParticle>,   p_cbCharged)
        ADD_BRANCH_T(int,           matchedsize)
        ADD_BRANCH_T(double,           angle_truerecon_eP)
        ADD_BRANCH_T(double,           angle_truerecon_eM)
        ADD_BRANCH_T(double,           angle_truerecon_Proton)
        ADD_BRANCH_T(TSimpleParticle,  recon_eM)
        ADD_BRANCH_T(TSimpleParticle,  recon_eP)
        ADD_BRANCH_T(TSimpleParticle,  recon_Proton)
        ADD_BRANCH_T(TSimpleParticle,  MCtrue_eM)
        ADD_BRANCH_T(TSimpleParticle,  MCtrue_eP)
        ADD_BRANCH_T(TSimpleParticle,  MCtrue_Proton)

        void fillAndReset() // to fill the tree and reset all values which add up
        {
            Tree->Fill();
            nTAPSneutral    = -1;
            nTAPScharged    = -1;
            nCBneutral      = -1;
            nCBcharged      = -1;
            nCombsInIM      = 0;
            IsMC            = false;
            nClusters       = 0;
            TaggW           = 0;
            p_tapsCharged().resize(0);
            p_cbCharged().resize(0);
            CB_effectiveradius().resize(0);
            CB_lateralmoment().resize(0);
            matchedsize = -1;
            angle_truerecon_eP = 0;
            angle_truerecon_eM = 0;
            angle_truerecon_Proton = 0;
//            recon_eM = TSimpleParticle();
//            recon_eP = TSimpleParticle();
//            recon_Proton = TSimpleParticle();
//            MCtrue_eM = TSimpleParticle();
//            MCtrue_eP = TSimpleParticle();
//            MCtrue_Proton = TSimpleParticle();
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
    TH1D* h_nCombsInIM = nullptr;

    TH2D* h_nCandCharged = nullptr;
    TH2D* h_nCand = nullptr;
    TH2D* h_cbdEE = nullptr;

    TH2D* h_clusteranalysis = nullptr;

    tree_t t;
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    double  CBAvgTime  = 0.0;
    double  CBSumVetoE  = 0.0;
    int hasnttriggered = 0;


private:

    utils::ClusterTools clustertools;

public:
    double effective_radius(const TCandidatePtr&) const;
    double lat_moment(const TCandidatePtr&) const;

};


}}} // end of ant::analysis::physics
