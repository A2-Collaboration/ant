#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/fitter/KinFitter.h"
#include "utils/uncertainties/Interpolated.h"
#include "utils/TriggerSimulation.h"

namespace ant {
namespace analysis {
namespace physics {

class FindCBESumThreshold : public Physics {
protected:

    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    TH1D* steps;


    TH1D* h_CBESum_raw;
    TH1D* h_CBESum_pr;

    TH1D* h_TaggT;
    TH1D* h_zvertex;



    using model_t = std::shared_ptr<const utils::UncertaintyModels::Interpolated>;

    model_t fit_model;
    utils::KinFitter fitter;


public:
    FindCBESumThreshold(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}}