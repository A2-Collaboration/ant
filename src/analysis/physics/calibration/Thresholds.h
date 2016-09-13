#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class Thresholds : public Physics {

protected:
    TH2D* hThresholds_ADC;
    TH2D* hThresholds_TDC;
    std::shared_ptr<Detector_t> Detector;

public:

    Thresholds(const Detector_t::Type_t& detectorType,
               const BinSettings& bins_x,
               const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics