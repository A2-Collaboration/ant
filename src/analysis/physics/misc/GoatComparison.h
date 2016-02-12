#pragma once

#include "analysis/physics/Physics.h"
#include <memory>

class TCanvas;

namespace ant {

class TH2CB;
class TH2TAPS;

namespace analysis {
namespace physics {

class GoatComparison : public Physics {
protected:

    const bool writeEvents;

    TH1D* steps;

    TH1D* h_CBSumVetoE;
    TH1D* h_PIDSumE;


    TH2CB* n_photon_high;
    TH2CB* n_photon_low;

    TH1D* IM_gg;

    TH2D* photon_thetaE;
    TH2D* photon_clusterSizeE;

public:
    GoatComparison(const std::string& name, OptionsPtr opts);
    virtual ~GoatComparison();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

};

}
}
}
