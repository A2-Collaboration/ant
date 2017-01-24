// better than using the standard ifndef include guards...
#pragma once

// physics classes need to derive from this interface
#include "physics/Physics.h"

// Ant provides many utility classes,
// such as support for prompt-random handling
#include "plot/PromptRandomHist.h"

// the physics classes reside in this nested namespace
namespace ant {
namespace analysis {
namespace physics {

// choose a nice name for your class
class Tutorial : public Physics {

    TH1D* h_nClusters;
    TH1D* h_nClusters_pr;


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