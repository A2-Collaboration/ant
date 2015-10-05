#pragma once

#include "analysis/physics/Physics.h"
class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {
namespace physics {

class EtapOmegaG : public Physics {

    BinSettings bins_im = BinSettings(1200);

    TH1D* gggg;
    TH1D* ggg;
    TH1D* gg;

    TH1D* IM_etap;
    TH1D* IM_omega;
    TH1D* IM_pi0;

    TH1D* Chi2_All;
    TH1D* Chi2_Min;

public:
    EtapOmegaG(PhysOptPtr opts);
    virtual void ProcessEvent(const data::Event& event) override;
    virtual void ShowResult() override;
};


}}}