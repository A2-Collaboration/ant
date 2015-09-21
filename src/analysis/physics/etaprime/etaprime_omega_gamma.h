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

public:
    EtapOmegaG(PhysOptPtr opts);
    virtual void ProcessEvent(const data::Event& event) override;
    virtual void ShowResult() override;
};


}}}