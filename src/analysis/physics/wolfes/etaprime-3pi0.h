#pragma once

#include "analysis/physics/Physics.h"
class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class Etap3pi0 : public Physics {

    TH1D* hNgamma;
    TH1D* hNgammaMC;

    TH1D* h2g;
    TH1D* h6g;

    TH1D* IM_etap;
    TH1D* IM_pi0;

    TH2D* IM_vs_chi2;

public:
    Etap3pi0(PhysOptPtr opts);
    virtual void ProcessEvent(const data::Event& event) override;
    virtual void ShowResult() override;
};


}}}
