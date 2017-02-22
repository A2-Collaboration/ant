#pragma once

#include "physics/Physics.h"
#include "utils/fitter/TreeFitter.h"
#include "base/WrapTTree.h"

namespace ant {
namespace analysis {
namespace physics {

class TestSigmaPlusFitter : public Physics {

    utils::TreeFitter treefitter;
    utils::TreeFitter::tree_t treefitter_SigmaPlus;
    utils::TreeFitter::tree_t treefitter_K0s;


    struct tree_t : WrapTTree {
        ADD_BRANCH_T(double, SigmaPlus_DeltaE)
        ADD_BRANCH_T(double, SigmaPlus_DeltaAngle)
        ADD_BRANCH_T(double, K0s_DeltaE)
        ADD_BRANCH_T(double, K0s_DeltaAngle)
    };

    tree_t t;

public:
    TestSigmaPlusFitter(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}} // end of ant::analysis::physics
