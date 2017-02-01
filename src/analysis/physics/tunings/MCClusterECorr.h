#pragma once

#include "physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class MCClusterECorr : public Physics {

    TH1D* h_Steps = nullptr;
    TH2D* h_LostMCTrue = nullptr;

    TH2D* h_nFills_CB = nullptr;
    TH2D* h_EtrueErec_CB = nullptr;

    TH2D* h_nFills_TAPS = nullptr;
    TH2D* h_EtrueErec_TAPS = nullptr;

public:
    MCClusterECorr(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};


}}} // end of ant::analysis::physics