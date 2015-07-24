#pragma once

#include "Calibration.h"
#include "expconfig/Detector_t.h"
#include "tree/TDataRecord.h"

#include <memory>

class TH1;

namespace ant {
namespace calibration {

class Timing : public Calibration::Module {

public:

    Timing(
            Detector_t::Type_t DetectorType,
            Calibration::Converter::ptr_t converter,
            const double defaultGain = 1.0, // default gain is 1.0
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

    const double DefaultGain;
    std::vector<double> Gains; // externally given, usually constant over one beamtime

    std::vector<double> Offsets; // to be calibrated values

};

}}  // namespace ant::calibration
