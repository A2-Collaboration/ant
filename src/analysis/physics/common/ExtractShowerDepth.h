#include "analysis/physics/Physics.h"

#include "base/WrapTTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class ExtractShowerDepth : public Physics {
public:

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(double, TrueZVertex)
        ADD_BRANCH_T(double, TrueTheta)
        ADD_BRANCH_T(double, TrueEk)
        ADD_BRANCH_T(double, Theta)
        ADD_BRANCH_T(double, Ek)
        ADD_BRANCH_T(double, ThetaCorr)
        ADD_BRANCH_T(double, ShowerDepth)
        ADD_BRANCH_T(double, RadiationLength)
    };

protected:

    const double param1;

    TH1D* steps;
    Tree_t t;


public:
    ExtractShowerDepth(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}}
