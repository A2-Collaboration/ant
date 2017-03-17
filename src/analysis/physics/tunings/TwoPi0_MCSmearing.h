#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/fitter/TreeFitter.h"
#include "base/ParticleTypeTree.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "base/WrapTTree.h"
#include "analysis/utils/PullsWriter.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class TwoPi0_MCSmearing : public Physics {

    TH1D* steps;
    PromptRandom::Switch promptrandom;

    TH2D* cb_pi0_channel    = nullptr;
    TH3D* cb_pi0_thetaE     = nullptr;
    TH3D* cb_pi0_EElement   = nullptr;
    TH2D* taps_pi0_channel  = nullptr;
    TH3D* taps_pi0_thetaE   = nullptr;
    TH3D* taps_pi0_EElement = nullptr;
    TH2D* cb_channel_correlation = nullptr;
    TH2D* taps_channel_correlation = nullptr;

    void FillIM(const TCandidate& c, double IM);

    void FillCorrelation(const TCandidate& c1, const TCandidate& c2);

public:
    TwoPi0_MCSmearing(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
