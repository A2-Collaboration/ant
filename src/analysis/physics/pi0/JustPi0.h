#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class JustPi0 : public Physics {
protected:
    TH1D* steps;


    PromptRandom::Switch promptrandom;
    PromptRandom::Hist1  h_missingmass;
    PromptRandom::Hist1  IM_2g;



public:
    JustPi0(const std::string& name, PhysOptPtr opts);

    virtual void ProcessEvent(const data::Event& event) override;
    virtual void ShowResult() override;
};

}}}
