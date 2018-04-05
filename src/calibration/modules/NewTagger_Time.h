#pragma once

#include "Time.h"
#include "expconfig/detectors/Tagger.h"

#include <memory>

namespace ant {

namespace calibration {

class NewTagger_Time : public Time
{

public:
    NewTagger_Time(std::shared_ptr<expconfig::detector::Tagger> tagg,
             std::shared_ptr<DataManager> calmgr,
             std::map<expconfig::detector::Tagger::TDCSector_t, Calibration::Converter::ptr_t> converters,
             double defaultOffset,
             std::shared_ptr<gui::PeakingFitFunction> fitFunction,
             const interval<double>& timeWindow
             );
};

}}  // namespace ant::calibration
