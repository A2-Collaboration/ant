#pragma once

#include "physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "cross_section.h"
#include "det_eff.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief The ppi0_2gamma class
 * Class uses basic kinematic cuts to find gp->pi0p (pi0->2g)
 * 
 */
class ppi0_2gamma: public Physics {
protected:
    struct steps_t : WrapTTree {

        ADD_BRANCH_T(bool, isSignal)
        ADD_BRANCH_T(std::string, cut)
        ADD_BRANCH_T(double, promptrandom)

        void AddStep(bool _isSignal, const std::string& _cut, double fillweight = 1)
        {
            isSignal = _isSignal;
            cut = _cut;
            promptrandom = fillweight;
            Tree->Fill();
        }
    };

    steps_t steps;
    utils::CrossSection cross_section;
    utils::DetEff detection_efficiency;
    PromptRandom::Switch promptrandom;

protected:
//    utils::KinFitter myfitter;
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

public:
    ppi0_2gamma(const std::string& name, OptionsPtr opts);
    virtual ~ppi0_2gamma() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};

}
}
}
