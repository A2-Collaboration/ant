#pragma once

#include "analysis/physics/Physics.h"
#include "utils/particle_tools.h"
class TH1D;
class TTree;

class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class McTrue3Pi0 : public Physics {

protected:
    ant::analysis::utils::ParticleVars proton;
    std::vector<ant::analysis::utils::ParticleVars> pi0s;

    TTree* mcTrue;

public:
    McTrue3Pi0(const std::string& name, PhysOptPtr opts);
    virtual void ProcessEvent(const data::Event& event) override;
    virtual void ShowResult() override;
};


}}}
