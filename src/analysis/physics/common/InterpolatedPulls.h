#include "analysis/physics/Physics.h"

#include "plot/PromptRandomHist.h"
#include "utils/PullsWriter.h"
#include "utils/Fitter.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class InterpolatedPulls : public Physics {
protected:

    PromptRandom::Switch promptrandom;

    TH1D* steps;

    utils::UncertaintyModelPtr model;
    utils::KinFitter fitter;

    utils::PullsWriter pullswriter;

public:
    InterpolatedPulls(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;
};


}}}
