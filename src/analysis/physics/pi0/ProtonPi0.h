#pragma once

#include "analysis/physics/Physics.h"

#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/fitter/KinFitter.h"
#include "base/WrapTTree.h"
#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace physics {

class ProtonPi0 : public Physics {
public:
    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(double, TaggW)
        ADD_BRANCH_T(double, TaggCh)
        ADD_BRANCH_T(double, TaggT)
        ADD_BRANCH_T(double, TaggE)
        ADD_BRANCH_T(double, CBAvgTime)

        ADD_BRANCH_T(double, FitProb)

        ADD_BRANCH_T(double, IM_2g)
        ADD_BRANCH_T(double, IM_2g_fitted)
        ADD_BRANCH_T(unsigned, nPhotonsCB)

        ADD_BRANCH_T(bool,   Proton_inCB)
        ADD_BRANCH_T(double, Proton_Ek)
        ADD_BRANCH_T(double, Proton_Theta)
        ADD_BRANCH_T(double, Proton_Phi)
        ADD_BRANCH_T(double, Proton_VetoE)
        ADD_BRANCH_T(double,   Proton_MinPIDPhi)
        ADD_BRANCH_T(unsigned, Proton_MinPIDCh)


        ADD_BRANCH_T(std::vector<unsigned>, PID_Ch)
        ADD_BRANCH_T(std::vector<double>, PID_Phi)
        ADD_BRANCH_T(std::vector<double>, PID_E)
        ADD_BRANCH_T(std::vector<double>, PID_Time)
    };

private:
    PromptRandom::Switch promptrandom;

    utils::KinFitter fitter;

    Tree_t t;

    TH1D* h_Steps;

public:
    ProtonPi0(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

};

}}}