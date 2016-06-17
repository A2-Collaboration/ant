#pragma once

#include "physics/Physics.h"
#include "utils/Fitter.h"
#include "plot/PromptRandomHist.h"

#include "base/ParticleTypeTree.h"

#include "root-addons/analysis_codes/hstack.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class FitPulls : public Physics {
protected:
    static const std::vector<ParticleTypeTreeDatabase::Channel> channels;

    struct ChannelHists_t {
        TH1D* h_probability;

        TH1D* p_cb_g_E = nullptr;
        TH1D* p_cb_g_Theta = nullptr;
        TH1D* p_cb_g_Phi = nullptr;
        TH1D* p_cb_p_Theta = nullptr;
        TH1D* p_cb_p_Phi = nullptr;
        TH1D* c_cb_g_Theta = nullptr;
        TH1D* c_cb_g_E = nullptr;
        TH1D* c_cb_p_Theta = nullptr;
        TH1D* p_taps_g_E = nullptr;
        TH1D* p_taps_g_Theta = nullptr;
        TH1D* p_taps_g_Phi = nullptr;
        TH1D* p_taps_p_Theta = nullptr;
        TH1D* p_taps_p_Phi = nullptr;
        TH1D* c_taps_g_Theta = nullptr;
        TH1D* c_taps_g_E = nullptr;
        TH1D* c_taps_p_Theta = nullptr;

        ChannelHists_t(const HistogramFactory& histFac);
    };

    struct ChannelItem_t {
        // use unique_ptr to make ctor easier...
        std::unique_ptr<utils::TreeFitter> Fitter;
        ParticleTypeTree                   Tree;
        std::unique_ptr<ChannelHists_t>    Hists_All;
        std::unique_ptr<ChannelHists_t>    Hists_Sig;
        unsigned          Multiplicity;
        ChannelItem_t(const HistogramFactory& parent,
                      ParticleTypeTreeDatabase::Channel channel,
                      utils::UncertaintyModelPtr uncertainty_model);
    };

    std::map<unsigned, std::vector<ChannelItem_t>> treefitters;

    bool findProton(const TCandidateList& cands, const TTaggerHit& taggerhit,
                    TParticlePtr& proton, TParticleList& photons,
                    LorentzVec& photon_sum);

    bool findProtonVetos(const TCandidateList& cands, const TTaggerHit& taggerhit,
                    TParticlePtr& proton, TParticleList& photons,
                    LorentzVec& photon_sum);

    PromptRandom::Switch promptrandom;

    TH1D* h_protoncopl;
    TH1D* h_taggtime;

    std::vector<ant::hstack*> hstacks_all;
    std::vector<ant::hstack*> hstacks_sig;


    const bool opt_save_after_cut;
    const bool opt_save_only;
public:
    FitPulls(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
