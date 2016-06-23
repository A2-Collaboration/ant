#pragma once

#include "Time.h"

#include <memory>

namespace ant {

namespace expconfig {
namespace detector {
struct TAPS;
}}

namespace calibration {

class TAPS_Time : public Time
{

public:
    TAPS_Time(std::shared_ptr<expconfig::detector::TAPS> taps,
              std::shared_ptr<DataManager> calmgr,
              Calibration::Converter::ptr_t converter_BaF2,
              Calibration::Converter::ptr_t converter_PbWO4,
              const interval<double>& timeWindow_BaF2 = {-std_ext::inf, std_ext::inf},
              const interval<double>& timeWindow_PbWO4 = {-std_ext::inf, std_ext::inf},
              std::shared_ptr<calibration::gui::PeakingFitFunction> fct = getDefaultFitFct()
              );

    static std::shared_ptr<calibration::gui::PeakingFitFunction> getDefaultFitFct();
};

}}  // namespace ant::calibration
