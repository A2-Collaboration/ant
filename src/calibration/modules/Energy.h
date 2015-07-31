#pragma once

#include "calibration/Calibration.h"
#include "expconfig/Detector_t.h"

#include "tree/TDataRecord.h" // for TKeyValue, TID

#include <memory>

class TH1;

namespace ant {

class CalibrationDataManager;

namespace calibration {

class Energy :
        public Calibration::PhysicsModule, // this makes this module abstract
        public ReconstructHook::DetectorReadHits
{

public:

    Energy(Detector_t::Type_t detectorType,
           std::shared_ptr<CalibrationDataManager> calmgr,
           Calibration::Converter::ptr_t converter,
           double defaultPedestal,
           double defaultGain,
           double defaultThreshold,
           double defaultRelativeGain);

    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits, extrahits_t& extrahits) override;

    // Updateable_traits interface
    virtual std::list<TID> GetChangePoints() const override { return {}; }
    void Update(const TID&) override {}

    virtual ~Energy();

protected:

    const Detector_t::Type_t DetectorType;

    std::shared_ptr<CalibrationDataManager> calibrationManager;

    const Calibration::Converter::ptr_t Converter;

    // only used for rawData conversion
    const double DefaultPedestal;
    std::vector<double> Pedestals;

    const double DefaultGain;
    std::vector<double> Gains;

    // always applied to values
    const double DefaultThreshold;
    std::vector<double> Thresholds;

    const double DefaultRelativeGain;
    std::vector<double> RelativeGains;

};

}}  // namespace ant::calibration
