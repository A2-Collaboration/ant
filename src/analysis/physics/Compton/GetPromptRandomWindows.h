#pragma once

#include "physics/Physics.h"
// To subtact out random tagger hits
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"

namespace ant {
namespace analysis {
namespace physics {

class GetPromptRandomWindows : public Physics {
public:
    GetPromptRandomWindows(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

    virtual void ShowResult() override;

private:
    TH1D* h_TaggerTime;
    TH1D* h_TaggerTimeCBSubtraction;

    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
};

}}}
