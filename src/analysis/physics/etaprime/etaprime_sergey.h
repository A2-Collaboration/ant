#include "analysis/physics/Physics.h"

#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/Fitter.h"
#include "analysis/utils/FitterSergey.h"
#include "base/WrapTTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class EtapSergey : public Physics {
public:

    struct params_t {
        const utils::UncertaintyModelPtr Fit_uncertainty_model;
        const bool Fit_Z_vertex;
        const double Z_vertex_sigma;
        params_t(utils::UncertaintyModelPtr fit_uncertainty_model,
                 bool fit_Z_vertex, double Z_vertex_sigma) :
            Fit_uncertainty_model(fit_uncertainty_model),
            Fit_Z_vertex(fit_Z_vertex),
            Z_vertex_sigma(Z_vertex_sigma)
        {}
    };

protected:

    const params_t params;

    utils::FitterSergey fitter_sergey;

    utils::KinFitter kinfitter_sig;

    utils::TreeFitter treefitter;
    utils::TreeFitter treefitter_Pi0Pi0;
    utils::TreeFitter treefitter_Pi0Eta;

    utils::TreeFitter::tree_t fitted_Pi0;
    utils::TreeFitter::tree_t fitted_g1_Pi0;
    utils::TreeFitter::tree_t fitted_g2_Pi0;

    utils::TreeFitter::tree_t fitted_Omega;
    utils::TreeFitter::tree_t fitted_g_Omega;

    utils::TreeFitter::tree_t fitted_EtaPrime;
    utils::TreeFitter::tree_t fitted_g_EtaPrime;

public:
    EtapSergey(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}}
