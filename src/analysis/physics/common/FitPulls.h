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
    static const std::vector<ParticleTypeTree> channels;

    std::map<unsigned, std::vector<utils::TreeFitter>> treefitters;

    void findProton(const TCandidateList& cands, const TTaggerHit& taggerhit,
                    TParticlePtr& proton, TParticleList& photons,
                    LorentzVec& photon_sum);

    PromptRandom::Switch promptrandom;

    TH1D* h_protoncopl;
    TH1D* h_taggtime;

    const bool opt_save_after_cut;
    const bool opt_save_only;
public:
    FitPulls(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}