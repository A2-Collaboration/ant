#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/Fitter.h"
#include "analysis/utils/FitterUncertainties.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/plot/PromptRandomHist.h"
#include "base/WrapTTree.h"

#include "root-addons/cbtaps_display/TH2CB.h"

#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace physics {

class PID_Energy_etaDalitz : public Physics {

protected:
    TH2* h_eegPID = nullptr;
    TH1D* h_counts = nullptr;
    TH1D* h_pTheta = nullptr;
    TH1D* h_protonVeto = nullptr;
    TH2* h_eta = nullptr;
    TH2* h_proton = nullptr;
    static constexpr double ETA_IM = 547.853;
    static constexpr double ETA_SIGMA = 50.;

    struct Tree_t : WrapTTree {
        Tree_t();

        ADD_BRANCH_T(TLorentzVector, eta)
        ADD_BRANCH_T(TLorentzVector, eta_fit)
        ADD_BRANCH_T(TLorentzVector, missing_momentum)
        ADD_BRANCH_T(double, copl)
        ADD_BRANCH_T(double, copl_final)

        ADD_BRANCH_T(double, TaggW)
        ADD_BRANCH_T(double, TaggE)
        ADD_BRANCH_T(double, TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)

        ADD_BRANCH_T(double, chi2)
        ADD_BRANCH_T(double, probability)
        ADD_BRANCH_T(unsigned, iterations)

        ADD_BRANCH_T(unsigned, nCands)
        ADD_BRANCH_T(double, CBSumE)
        ADD_BRANCH_T(double, CBAvgTime)
    };

    struct PerChannel_t {
        std::string title;
        TH2* eegPID = nullptr;
        TH1D* steps = nullptr;
        TH1D* etaIM = nullptr;
        TH1D* etaIM_fit = nullptr;
        TH1D* etaIM_cand = nullptr;
        TH1D* etaIM_final = nullptr;
        TH1D* MM = nullptr;
        TH1D* hCopl = nullptr;
        TH1D* hCopl_final = nullptr;
        TH1D* hChi2 = nullptr;
        TH1D* hProb = nullptr;
        TH1D* hIter = nullptr;

        TH2* proton_E_theta = nullptr;

        PerChannel_t(const std::string& Title, HistogramFactory& hf);

        void Show();
        void Fill(const TEventData& d);
    };

    std::map<std::string, PerChannel_t> channels;
    Tree_t t;
    PromptRandom::Switch promptrandom;
    utils::KinFitter kinfit;
    using uncertainty_model_t = utils::UncertaintyModels::Optimized_Oli1;

    template<typename T>
    void shift_right(std::vector<T>&);

public:

    PID_Energy_etaDalitz(const std::string& name, OptionsPtr opts);

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
