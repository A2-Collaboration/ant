#pragma once

#include "analysis/physics/Physics.h"
#include "base/piecewise_interval.h"
#include "base/WrapTTree.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {
namespace physics {

class ReconstructCheck : public Physics {

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(unsigned,            Multiplicity)
        ADD_BRANCH_T(std::vector<double>, Thetas)
        ADD_BRANCH_T(std::vector<double>, Energies)
    };

    Tree_t t;

    TH1D* hMultiplicities = nullptr;

    PiecewiseInterval<unsigned> Multiplicities;

public:
    ReconstructCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
