#pragma once

#include "analysis/physics/Physics.h"

#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/fitter/KinFitter.h"

namespace ant {
namespace analysis {
namespace physics {

class ThreePhotonCheck : public Physics {

    PromptRandom::Switch promptrandom;

    utils::KinFitter fitter;

    TH1D* h_Steps;
    TH1D* h_photonsIM;

public:
    ThreePhotonCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

};

}}}
