#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/Fitter.h"
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
public:

public:
    TwoPi0_MCSmearing(const std::string& name, OptionsPtr opts);
    virtual ~TwoPi0_MCSmearing();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

};

}}}
