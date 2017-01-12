#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class PID_Energy_Banana : public Physics {

protected:
    TH2D* h_pedestals = nullptr;
    TH3D* h_bananas = nullptr;

    struct PerChannel_t {
        TH2D* PedestalTiming = nullptr;
        TH1D* PedestalNoTiming = nullptr;
        TH2D* Banana = nullptr;
        TH2D* BananaRaw = nullptr;
        TH1D* TDCMultiplicity;
        TH1D* QDCMultiplicity;
        //TH3D* BananaTiming = nullptr;
        PerChannel_t(analysis::HistogramFactory HistFac);
    };

    std::vector<PerChannel_t> h_perChannel;
public:

    PID_Energy_Banana(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
