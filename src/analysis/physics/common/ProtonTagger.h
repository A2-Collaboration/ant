#pragma once

#include "analysis/physics/Physics.h"
#include "base/interval.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {
namespace physics {

class ProtonTagger : public Physics {
public:

    //const interval<double> pi0_cut = interval<double>::CenterWidth(ParticleTypeDatabase::Pi0.Mass(),    10.0);
    const interval<double> pi0_cut   = interval<double>::CenterWidth(126, 10);
    const interval<double> mm_cut  = interval<double>::CenterWidth(1025, 50.0);

    TH2D* tof = nullptr;
    TH2D* dEE = nullptr;
    TH2D* cls = nullptr;

    TH1D* ggIM = nullptr;
    TH1D* MM_after_cut = nullptr;
    TH1D* angle = nullptr;


    ProtonTagger(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

}
}
}
