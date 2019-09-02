#pragma once

#include "analysis/physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "analysis/utils/TriggerSimulation.h"
#include "analysis/plot/PromptRandomHist.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A class for analyzing the pi0->e+e-g decay
 *
 */
class Pi0Dalitz: public Physics {
public:
    struct AnalysTree_t : WrapTTree {

        //--- Trigger
        ADD_BRANCH_T(bool, TBHasTrigg)
        //--- MCstuff
        ADD_BRANCH_T(bool, TBMCpi0eeg)
        ADD_BRANCH_T(TLorentzVector, TBTrueVecpi0)
        ADD_BRANCH_T(TLorentzVector, TBTrueVecep)
        ADD_BRANCH_T(TLorentzVector, TBTrueVecem)
        ADD_BRANCH_T(std::vector<TLorentzVector>, TBTrueVecGammas)
        //--- All candidate info
        ADD_BRANCH_T(std::vector<double>, TBCanCaloE)
        ADD_BRANCH_T(std::vector<double>, TBCanTheta)
        ADD_BRANCH_T(std::vector<double>, TBCanPhi)
        ADD_BRANCH_T(std::vector<double>, TBCanTime)
        ADD_BRANCH_T(std::vector<double>, TBCanCluSize)
        ADD_BRANCH_T(std::vector<bool>, TBCanInCB)
        //--- Info from tagger
        ADD_BRANCH_T(std::vector<double>, TBPrRndWeight)
        ADD_BRANCH_T(std::vector<TLorentzVector>, TBInitPhotVec)
    };


private:
    AnalysTree_t AnalysTree;
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    void CreateHistos();

public:
    Pi0Dalitz(const std::string& name, OptionsPtr opts);
    virtual ~Pi0Dalitz() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};

}
}
}
