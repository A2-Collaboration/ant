#pragma once

#include "analysis/physics/Physics.h"
#include "utils/ParticleTools.h"
class TH1D;
class TTree;

class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class McTrue3Pi0 : public Physics {

protected:

    const std::vector<std::vector<unsigned>> combinations =
    {
        { 0, 1, 2, 3, 4, 5 },
        { 0, 1, 2, 4, 3, 5 },
        { 0, 1, 2, 5, 3, 4 },

        { 0, 2, 1, 3, 4, 5 },
        { 0, 2, 1, 4, 3, 5 },
        { 0, 2, 1, 5, 3, 4 },

        { 0, 3, 1, 2, 4, 5 },
        { 0, 3, 1, 4, 2, 5 },
        { 0, 3, 1, 5, 2, 4 },

        { 0, 4, 1, 2, 3, 5 },
        { 0, 4, 1, 3, 2, 5 },
        { 0, 4, 1, 5, 2, 3 },

        { 0, 5, 1, 2, 3, 4 },
        { 0, 5, 1, 3, 2, 4 },
        { 0, 5, 1, 4, 2, 3 }
    };

    ant::analysis::utils::ParticleVars proton;
    std::vector<ant::analysis::utils::ParticleVars> pi0s;
    std::vector<double> popens;
    std::vector<ant::analysis::utils::ParticleVars> gammas;

    std::vector<double> GetAllPhotonAngles(const TParticleList& photons) const;

    TH1D* hAngle;
    TTree* mcTrue;

public:
    McTrue3Pi0(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}}
