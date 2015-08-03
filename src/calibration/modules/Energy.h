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
        public Calibration::Module, // this makes this module abstract
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
    virtual std::vector<std::list<TID>> GetChangePoints() const override;
    void Update(std::size_t index, const TID& tid) override;

    virtual ~Energy();

protected:

    const Detector_t::Type_t DetectorType;

    std::shared_ptr<CalibrationDataManager> calibrationManager;

    const Calibration::Converter::ptr_t Converter;

    struct CalibType
    {
        const double        DefaultValue;
        std::vector<double> Values;
        const std::string   Type;
        const std::size_t   Index;
        CalibType(double defaultValue, const std::string type, size_t index):
            DefaultValue(defaultValue),
            Values(),
            Type(type),
            Index(index)
        {}
    };


    CalibType Pedestals;
    CalibType Gains;
    CalibType Thresholds;
    CalibType RelativeGains;

};

}}  // namespace ant::calibration
