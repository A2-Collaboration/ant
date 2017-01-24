// better than using the standard ifndef include guards...
#pragma once

// physics classes need to derive from this interface
#include "physics/Physics.h"

// Ant provides many utility classes,
// such as support for prompt-random handling
#include "plot/PromptRandomHist.h"

// handle ROOT's TTree
#include "base/WrapTTree.h"

// the physics classes reside in this nested namespace
namespace ant {
namespace analysis {
namespace physics {

// choose a nice name for your class
class Tutorial : public Physics {
public:

    // use "struct" to avoid all those public keywords...
    struct tree_t : WrapTTree {
        // define two branches with ADD_BRANCH_T (again, this is black macro magic)
        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(unsigned, nClusters)
        // you may even use more complex types with non-default ctor, for example
//        ADD_BRANCH_T(std::vector<double>, Numbers, 3) // vector Numbers is three items large now
    };


private:

    TH1D* h_nClusters;
    TH1D* h_nClusters_pr;

    tree_t t;

    // use this instance to handle prompt-random weighting
    PromptRandom::Switch promptrandom;

public:
    // physics need to implement this public constructor...
    Tutorial(const std::string& name, OptionsPtr opts);

    // ...and the following method (use override keyword)
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

    // optional method (called in interactive mode to show histograms)
    virtual void ShowResult() override;
};


}}} // end of ant::analysis::physics