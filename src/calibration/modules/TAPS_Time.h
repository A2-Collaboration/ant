#pragma once

#include "Time.h"


namespace ant {

namespace expconfig {
namespace detector {
class TAPS;
}}

namespace calibration {

class TAPS_Time : public Time
{

public:
    TAPS_Time(std::shared_ptr<expconfig::detector::TAPS> taps,
              std::shared_ptr<DataManager> calmgr,
              Calibration::Converter::ptr_t converter,
              const interval<double>& timeWindow = {-std_ext::inf, std_ext::inf}
              );
};

}}  // namespace ant::calibration
