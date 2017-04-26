#include "analysis/physics/Physics.h"

#include "plot/PromptRandomHist.h"
#include "utils/PullsWriter.h"
#include "utils/fitter/KinFitter.h"
#include "utils/uncertainties/Interpolated.h"
#include "analysis/utils/TriggerSimulation.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {
namespace physics {

class InterpolatedPulls : public Physics {
protected:

    utils::TriggerSimulation triggersimu;
    PromptRandom::Switch promptrandom;

    TH1D* steps;
    TH1D* h_missingmass_cut;
    TH1D* h_missingmass_best;
    TH2D* h_IM_gg_gg;
    TH2D* h_IM_gg_gg_cut;

    TH1D* h_zvertex;

    TH2D* h_proton_E_theta;

    TH2D* h_ToF_E_photon_taps;
    TH2D* h_ToF_E_proton_taps;
    TH2D* h_PSA_photon_taps;
    TH2D* h_PSA_proton_taps;

    TH2D* h_E_vetoE_photon_cb;
    TH2D* h_E_vetoE_photon_taps;
    TH2D* h_E_vetoE_proton_cb;
    TH2D* h_E_vetoE_proton_taps;

    using model_t = std::shared_ptr<const utils::UncertaintyModels::Interpolated>;

    model_t fit_model;
    utils::KinFitter fitter;

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(std::vector<double>, IM_gg_gg, 2)
        ADD_BRANCH_T(double, IM_pi0pi0_radius)
    };

    utils::PullsWriter<Tree_t> pullswriter;

    Tree_t t;

public:
    InterpolatedPulls(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;
};


}}}
