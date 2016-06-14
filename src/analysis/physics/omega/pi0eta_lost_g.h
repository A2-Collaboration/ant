#pragma once
#include "analysis/physics/Physics.h"
#include "base/Tree.h"
#include "base/interval.h"
#include "base/std_ext/math.h"
#include "base/WrapTTree.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Rtypes.h"
#pragma GCC diagnostic pop

#include <map>


class TH1D;
class TH2D;
class TH3D;

namespace ant {

namespace analysis {
namespace physics {



class Pi0EtaLostG : public Physics {
public:
    struct Tree_t : WrapTTree {
        Tree_t();

//        ADD_BRANCH_T(std::vector<TLorentzVector>, photon, 3)
//        ADD_BRANCH_T(TLorentzVector, proton)
        ADD_BRANCH_T(double, E)
        ADD_BRANCH_T(double, theta)
        ADD_BRANCH_T(double, phi)
        ADD_BRANCH_T(double, MinAngle)
    };

protected:
    TTree*  tree;
    Tree_t t;

public:

    Pi0EtaLostG(const std::string& name, OptionsPtr opts);
    virtual ~Pi0EtaLostG();

    void ProcessEvent(const TEvent& event, manager_t&) override;
};

}
}
}
