#pragma once

#include "analysis/physics/Physics.h"
class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {
namespace physics {

class EtapOmegaG : public Physics {

    struct perDecayHists_t {
        TH1D* gggg;
        TH1D* ggg;
        TH1D* gg;

        TH2D* IM_etap_omega;
        TH1D* IM_pi0;

        TH1D* Chi2_All;
        TH1D* Chi2_Best;

        perDecayHists_t(SmartHistFactory& HistFac, const std::string& decaystring);
    };

    std::map<std::string, perDecayHists_t> perDecayHists;


public:
    EtapOmegaG(PhysOptPtr opts);
    virtual void ProcessEvent(const data::Event& event) override;
    virtual void ShowResult() override;
};


}}}