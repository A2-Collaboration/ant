#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/Fitter.h"
#include "utils/MCSmear.h"
#include "base/ParticleTypeTree.h"
#include "TLorentzVector.h"

class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class KinFitPi0 : public Physics {
protected:

    struct MultiPi0 {
        MultiPi0(HistogramFactory& histFac, unsigned nPi0, std::shared_ptr<const utils::Fitter::UncertaintyModel> fitter_model=utils::UncertaintyModels::MCExtracted::makeAndLoad());

        void ProcessData(const TEventData& data, const utils::MCSmear& smear);
        void ShowResult();

    protected:
        const unsigned multiplicity;

        HistogramFactory HistFac;


        using IM_perms_t = std::vector<std::vector<size_t>>;
        const IM_perms_t IM_perms;
        static IM_perms_t BuildIMPerms(unsigned multiplicity);

        utils::KinFitter fitter;

        TH1D* steps;

        TH1D* Proton_Coplanarity;
        TH1D* h_taggtime;

        TTree* tree = nullptr;
        TLorentzVector b_g1;
        TLorentzVector b_g2;
        TLorentzVector b_p;
        double         b_tagw;

        PromptRandom::Switch promptrandom;

        PromptRandom::Hist1 h_missingmass;
        PromptRandom::Hist1 h_fitprobability;
        PromptRandom::Hist1 IM_2g_byFit;
        PromptRandom::Hist1 IM_2g_fitted;

        std::map<std::string, PromptRandom::Hist1> h_pulls;
        using fit_variables_t = decltype(APLCON::Result_t::Variables);
        void FillPulls(const fit_variables_t& vars);
    };

    std::vector<std::unique_ptr<MultiPi0>> multiPi0;

    std::shared_ptr<utils::Fitter::UncertaintyModel> getModel(const std::string& model_name) const;

    std::shared_ptr<utils::Fitter::UncertaintyModel> fitter_model;
    std::shared_ptr<utils::Fitter::UncertaintyModel> smear_model;
    utils::MCSmear smear;
    const bool opt_useMCSmear = false;

public:
    KinFitPi0(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
