#pragma once

#include "physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "utils/TriggerSimulation.h"
#include "analysis/plot/PromptRandomHist.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A gp->ppi0 class
 *
 */
class scratch_lheijken_gppi0p: public Physics {
protected:
    struct AnalysTree_t : WrapTTree {

        //--- Checks
        ADD_BRANCH_T(double, TBPrRndWeight)
        //--- MCstuff
        ADD_BRANCH_T(bool, TBMCpi0gg)
        ADD_BRANCH_T(std::vector<TLorentzVector>, TBTrueGammas)
        ADD_BRANCH_T(TLorentzVector, TBTrueVecpi0)

        //--- All candidate info
        ADD_BRANCH_T(std::vector<double>, TBNeuCanCaloE)
        ADD_BRANCH_T(std::vector<double>, TBNeuCanTheta)
        ADD_BRANCH_T(std::vector<double>, TBNeuCanPhi)
        ADD_BRANCH_T(std::vector<double>, TBNeuCanTime)
        ADD_BRANCH_T(std::vector<double>, TBNeuCanCluSize)
        ADD_BRANCH_T(std::vector<bool>, TBNeuCanInCB)
        ADD_BRANCH_T(std::vector<double>, TBChaCanCaloE)
        ADD_BRANCH_T(std::vector<double>, TBChaCanVetoE)
        ADD_BRANCH_T(std::vector<double>, TBChaCanTheta)
        ADD_BRANCH_T(std::vector<double>, TBChaCanPhi)
        ADD_BRANCH_T(std::vector<double>, TBChaCanTime)
        ADD_BRANCH_T(std::vector<double>, TBChaCanCluSize)
        ADD_BRANCH_T(std::vector<bool>, TBChaCanInCB)
        ADD_BRANCH_T(std::vector<bool>, TBChaCanInPID)
        //--- pi0gg stuff
        ADD_BRANCH_T(std::vector<double>, TBpi0ggTime)
        ADD_BRANCH_T(std::vector<bool>, TBpi0ggOnlyCB)
        ADD_BRANCH_T(std::vector<TLorentzVector>, TBpi0ggVecRec)
        ADD_BRANCH_T(std::vector<int>, TBpi0ggNeuCombInd1)
        ADD_BRANCH_T(std::vector<int>, TBpi0ggNeuCombInd2)
        //--- pi0gg stuff
        ADD_BRANCH_T(std::vector<double>, TBpi0DDTime)
        ADD_BRANCH_T(std::vector<bool>, TBpi0DDOnlyCB)
        ADD_BRANCH_T(std::vector<TLorentzVector>, TBpi0DDVecRec)


    };
    AnalysTree_t AnalysTree;
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    //--- MC True
    TH1D *hTrueGammaE, *hTrueIMgg, *hTrueMMgg;
    TH2D *hTrueThevsEg, *hTrueThevsPhig, *hTrueThevsEpi0, *hTrueThevsPhipi0;
    //--- Checks
    TH1D *hPromRandWei;
    //--- All candidate info
    TH2D *hNeuCanThevsCaloE, *hNeuCanThevsPhi;
    TH1D *hNeuCanCBTime, *hNeuCanTAPSTime, *hNeuCanCBCluSize, *hNeuCanTAPSCluSize;
    TH2D *hChaCanThevsCaloE,*hChaCanThevsVetoE, *hChaCanThevsPhi;
    TH1D *hChaCanCBTime, *hChaCanTAPSTime, *hChaCanCBCluSize, *hChaCanTAPSCluSize;
    //--- pi0gg stuff
    TH1D *hpi0ggO2gOCB_time,*hpi0ggO2gCBTA_time, *hpi0ggO2gOCB_IM, *hpi0ggO2gCBTA_IM;
    TH2D *hpi0ggO2g_gThevsE, *hpi0ggO2g_gThevsPhi, *hpi0ggO2g_pi0ThevsE, *hpi0ggO2g_pi0ThevsPhi;
    //--- pi0DD stuff
    TH1D *hpi0DD1N2COCB_time,*hpi0DD1N2CCBTA_time, *hpi0DD1N2COCB_IM, *hpi0DD1N2CCBTA_IM;
    TH2D *hpi0DD1N2C_pi0ThevsE, *hpi0DD1N2C_pi0ThevsPhi;

protected:
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    void CreateHistos();

public:
    scratch_lheijken_gppi0p(const std::string& name, OptionsPtr opts);
    virtual ~scratch_lheijken_gppi0p() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};

}
}
}
