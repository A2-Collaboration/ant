#pragma once

#include "calibration/Calibration.h"
#include "base/Detector_t.h"


namespace ant {
namespace calibration {

class Scaler :
        public Calibration::BaseModule,
        public ReconstructHook::DetectorReadHits
{

public:

    Scaler(Detector_t::Type_t DetectorType,
           Calibration::Converter::ptr_t converter);
    virtual ~Scaler();

    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits) override;

protected:
    const Detector_t::Type_t DetectorType;
    const Calibration::Converter::ptr_t Converter;
};

}}  // namespace ant::calibration
