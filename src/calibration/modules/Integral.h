#pragma once

#include "Calibration.h"
#include "expconfig/Detector_t.h"
#include "tree/TDataRecord.h"

#include <memory>

class TH1;

namespace ant {
namespace calibration {

class Integral : public Calibration::Module {

public:

    Integral(Detector_t::Type_t DetectorType,
             Calibration::Converter::ptr_t converter,
             const double defaultPedestal,
             const double defaultGain,
             const double defaultThreshold);

    // CalibrationApply_traits interface
    virtual void ApplyTo(const readhits_t& hits, extrahits_t& extrahits) override;
    virtual void ApplyTo(event_ptr&) override {}


    // Physics interface
    void ProcessEvent(const Event &event) override;
    void Finish() override;
    void ShowResult() override;

    // CalibrationUpdate_traits interface
    void BuildRanges(std::list<TID>&) override {}
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
