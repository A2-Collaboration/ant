#pragma once

#include "Calibration.h"
#include "expconfig/Detector_t.h"
#include "tree/TDataRecord.h"
#include "base/interval.h"

#include <memory>
#include <limits>

class TH1;

namespace ant {
namespace calibration {

class Timing : public Calibration::Module {

public:

    Timing(
            Detector_t::Type_t DetectorType,
            Calibration::Converter::ptr_t converter,
            double defaultOffset,
            const interval<double>& timeWindow = {-std_ext::inf, std_ext::inf},
            double defaultGain = 1.0, // default gain is 1.0
            const std::vector< TKeyValue<double> >& gains = {}
            );

    // CalibrationApply_traits interface
    virtual void ApplyTo(const std::map< Detector_t::Type_t, std::list< TDetectorReadHit* > >& hits) override;
    virtual void ApplyTo(std::unique_ptr<TEvent>&) override {}


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

    const interval<double> TimeWindow;

    const double DefaultOffset;
    std::vector<double> Offsets;

    const double DefaultGain;
    std::vector<double> Gains;

};

}}  // namespace ant::calibration
