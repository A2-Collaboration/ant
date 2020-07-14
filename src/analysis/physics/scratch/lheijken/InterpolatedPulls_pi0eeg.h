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

class scratch_lheijken_InterpolatedPulls_pi0eeg : public Physics {
protected:

    utils::TriggerSimulation triggersimu;
    PromptRandom::Switch promptrandom;

    TH1D* steps;
    TH1D* h_missingmass_cut;
    TH1D* h_missingmass_best;
    TH1D* h_IM_eeg;
    TH1D* h_IM_eeg_cut;

    TH1D* h_zvertex;

    TH2D* h_proton_E_theta;

    TH2D* h_E_vetoE_photon_cb;
    TH2D* h_E_vetoE_photon_taps;
    TH2D* h_E_vetoE_proton_cb;
    TH2D* h_E_vetoE_proton_taps;

    using model_t = std::shared_ptr<const utils::UncertaintyModels::Interpolated>;

    model_t fit_model_data;
    model_t fit_model_mc;
    utils::KinFitter fitter;

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(double, IM_egg)
        ADD_BRANCH_T(double, IM_pi0_radius)
    };

    utils::PullsWriter<Tree_t> pullswriter;

    Tree_t t;

public:
    scratch_lheijken_InterpolatedPulls_pi0eeg(const std::string& name, OptionsPtr opts);
    virtual ~scratch_lheijken_InterpolatedPulls_pi0eeg() {}
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;
};


}}}
