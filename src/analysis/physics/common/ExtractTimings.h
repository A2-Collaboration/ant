#pragma once

#include "analysis/physics/Physics.h"
#include "base/interval.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {
namespace physics {

class ExtractTimings : public Physics {
public:

    TH1D* EPT_CB   = nullptr;
    TH1D* EPT_TAPS = nullptr;
    TH1D* CB_TAPS  = nullptr;
    TH1D* CB_CB    = nullptr;
    TH1D* TAPS_TAPS  = nullptr;
    TH1D* CB       = nullptr;
    TH1D* TAPS     = nullptr;
    TH1D* EPT      = nullptr;

    ExtractTimings(const std::string& name, PhysOptPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}
}
}
