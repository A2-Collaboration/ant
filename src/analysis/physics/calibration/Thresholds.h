#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class Thresholds : public Physics {

protected:
    TH2D* hThresholds_Raw;
    TH2D* hThresholds_ADC;
    TH2D* hThresholds_TDC;
    TH2D* hThresholds_Raw_TDC;
    TH1D* hADCnoTDC;
    TH1D* hADC;
    TH1D* hADCnoTDC_norm = nullptr;
    const std::shared_ptr<const Detector_t> Detector;

    const double badTDCthreshold;

public:

    Thresholds(const Detector_t::Type_t& detectorType,
               const BinSettings& bins_x,
               const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;
};

}}} // namespace ant::analysis::physics
