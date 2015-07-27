#pragma once

#include "Calibration.h"
#include "expconfig/Detector_t.h"

#include "tree/TDataRecord.h" // for TKeyValue, TID

#include <memory>

class TH1;

namespace ant {
namespace calibration {

class Integral : public Calibration::Module {

public:

    Integral(Detector_t::Type_t detectorType,
             Calibration::Converter::ptr_t converter,
             const double defaultPedestal,
             const double defaultGain,
             const double defaultThreshold);

    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits, extrahits_t& extrahits) override;

    // Physics interface
    void ProcessEvent(const Event &event) override;
    void Finish() override;
    void ShowResult() override;

    // Updateable_traits interface
    virtual std::list<TID> GetChangePoints() const override { return {}; }
    void Update(const TID&) override {}

protected:

    const Detector_t::Type_t DetectorType;
    const Calibration::Converter::ptr_t Converter;

    const double DefaultPedestal;
    std::vector<double> Pedestals;

    const double DefaultGain;
    std::vector<double> Gains;

    const double DefaultThreshold;
    std::vector<double> Thresholds;

};

}}  // namespace ant::calibration
