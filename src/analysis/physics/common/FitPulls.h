#pragma once

#include "physics/Physics.h"
#include "utils/Fitter.h"
#include "plot/PromptRandomHist.h"

#include "base/ParticleTypeTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class FitPulls : public Physics {
protected:
    static const std::vector<ParticleTypeTreeDatabase::Channel> channels;

    struct ChannelHists_t {
        TH1D* p_cb_g_E;
        TH1D* p_cb_g_Theta;
        TH1D* p_cb_g_Phi;
        TH1D* p_cb_p_E;
        TH1D* p_cb_p_Theta;
        TH1D* p_cb_p_Phi;
        TH1D* p_taps_g_E;
        TH1D* p_taps_g_Theta;
        TH1D* p_taps_g_Phi;
        TH1D* p_taps_p_E;
        TH1D* p_taps_p_Theta;
        TH1D* p_taps_p_Phi;

        ChannelHists_t(const HistogramFactory& histFac, const std::string& name);
    };

    struct ChannelItem_t {
        // use unique_ptr to make ctor easier...
        std::unique_ptr<utils::TreeFitter> Fitter;
        std::unique_ptr<ChannelHists_t>    Hists;
        unsigned          Multiplicity;
        ChannelItem_t(const HistogramFactory& parent,
                      ParticleTypeTreeDatabase::Channel channel,
                      utils::Fitter::UncertaintyModelPtr uncertainty_model);
    };

    std::map<unsigned, std::vector<ChannelItem_t>> treefitters;

    void findProton(const TCandidateList& cands, const TTaggerHit& taggerhit,
                    TParticlePtr& proton, TParticleList& photons,
                    LorentzVec& photon_sum);

    PromptRandom::Switch promptrandom;

    TH1D* h_protoncopl;
    TH1D* h_taggtime;
    TH1D* h_probability;

    const bool opt_save_after_cut;
    const bool opt_save_only;
public:
    FitPulls(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
