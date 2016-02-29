#pragma once

#include "analysis/physics/Physics.h"

#include "root-addons/cbtaps_display/TH2CB.h"

namespace ant {
namespace analysis {
namespace physics {

class PID_Energy_etaDalitz : public Physics {

protected:
    TH2* eegPID = nullptr;
    static constexpr double ETA_IM = 547.853;
    static constexpr double ETA_SIGMA = 25.;

    template<typename T>
    void shift_right(std::vector<T>&);

public:

    PID_Energy_etaDalitz(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
