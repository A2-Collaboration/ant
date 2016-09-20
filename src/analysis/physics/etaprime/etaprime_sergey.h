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

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(double, TaggE)
        ADD_BRANCH_T(double, TaggT)
        ADD_BRANCH_T(double, TaggCh)

        ADD_BRANCH_T(double, KinFitProb)
        ADD_BRANCH_T(double, TreeFitProb)
        ADD_BRANCH_T(double, AntiPi0FitProb)
        ADD_BRANCH_T(double, AntiEtaFitProb)

        ADD_BRANCH_T(std::vector<double>, IM_3g)
        ADD_BRANCH_T(double, IM_4g)

        ADD_BRANCH_T(std::vector<double>, gNonPi0_Theta)
        ADD_BRANCH_T(std::vector<double>, gNonPi0_CaloE)

        ADD_BRANCH_T(double, CBVetoSumE)
        ADD_BRANCH_T(double, PIDSumE)

        ADD_BRANCH_T(unsigned, KinFitProtonIdx)
        ADD_BRANCH_T(unsigned, TreeFitProtonIdx)
        ADD_BRANCH_T(unsigned, AntiPi0FitProtonIdx)
        ADD_BRANCH_T(unsigned, AntiEtaFitProtonIdx)

        ADD_BRANCH_T(unsigned, MCTrue)

    };

    Tree_t treeSergey;
    Tree_t treeAnt;


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

    using result_t = utils::FitterSergey::result_t;

    const params_t params;

    utils::FitterSergey fitter_sergey;

    utils::KinFitter kinfitter;

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

    static void fillTree(Tree_t& t, const std::vector<result_t>& results,
                         unsigned MCTrue, double PIDSumE);

    TH1D* h_MissedBkg;
    static const ParticleTypeTree ptreeSignal;
    struct Background_t {
        const std::string Name;
        const ParticleTypeTree Tree;
        Background_t(const std::string& name, ParticleTypeTree tree) :
            Name(name), Tree(tree) {}
    };

public:
    EtapSergey(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

    static const std::vector<Background_t> ptreeBackgrounds;
};


}}}
