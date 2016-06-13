#pragma once
#include "analysis/physics/Physics.h"
#include "base/ParticleTypeTree.h"

class TH1D;
class TH2D;
class TH3D;

namespace ant {

namespace analysis {
namespace physics {

class OmegaMCAngles : public Physics {
public:
    OmegaMCAngles(const std::string& name, OptionsPtr opts);
    virtual ~OmegaMCAngles();

    void ProcessEvent(const TEvent& event, manager_t&) override;
    void Finish() override;
    void ShowResult() override;

protected:
    using decaytree_t = ant::Tree<const ParticleTypeDatabase::Type&>;

    const std::shared_ptr<decaytree_t> channelTree;


    TH1D* angle_pi0_gg  = nullptr;
    TH1D* angle_eta_gg  = nullptr;
    TH1D* angle_all     = nullptr;
    TH1D* angle_min_all = nullptr;
};

}
}
}
