#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/fitter/KinFitter.h"
#include "utils/Uncertainties.h"
#include "utils/TriggerSimulation.h"

namespace ant {
namespace analysis {
namespace physics {

class TriggerSimulation : public Physics {
protected:
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    TH1D* steps;

    TH1D* h_CBESum_raw;
    TH1D* h_CBESum_pr;

    TH1D* h_TaggT;
    TH1D* h_zvertex;


    utils::UncertaintyModelPtr fit_model;
    utils::KinFitter fitter;


public:
    TriggerSimulation(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}}