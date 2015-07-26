#pragma once

#include "Calibration.h"
#include "expconfig/Detector_t.h"


namespace ant {
namespace calibration {

class Scaler : public Calibration::SimpleModule {

public:

    Scaler(Detector_t::Type_t DetectorType,
           Calibration::Converter::ptr_t converter);
    virtual ~Scaler();

    // CalibrationApply_traits interface
    virtual void ApplyTo(const readhits_t& hits, extrahits_t&) override;
    virtual void ApplyTo(event_ptr&) override {}


protected:
    const Detector_t::Type_t DetectorType;
    const Calibration::Converter::ptr_t Converter;
};

}}  // namespace ant::calibration
