#pragma once

#include "calibration/Calibration.h"
#include "base/interval.h"

namespace ant {
namespace calibration {

class Tagger_QDC :
        public ReconstructHook::DetectorReadHits
{
public:
    Tagger_QDC(Detector_t::Type_t detectorType,
               Calibration::Converter::ptr_t converter);
    virtual ~Tagger_QDC();

    virtual void ApplyTo(const readhits_t& hits) override;
protected:
    const Detector_t::Type_t DetectorType;
    const Calibration::Converter::ptr_t Converter;
};

}}
