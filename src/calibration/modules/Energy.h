#pragma once

#include "CalibType.h"

#include "calibration/Calibration.h"
#include "base/Detector_t.h"

#include <memory>

class TH1;

namespace ant {

namespace calibration {

class Energy :
        public Calibration::Module, // this makes this module abstract
        public ReconstructHook::DetectorReadHits
{

public:
    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits) override;

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;
    void UpdatedTIDFlags(const TID& id) override;

protected:

    bool IsMC = false; // managed by UpdatedTIDFlags

    using defaults_t = const std::vector<double>&;

    Energy(const detector_ptr_t& det,
           const std::shared_ptr<DataManager>& calmgr,
           const Calibration::Converter::ptr_t& converter,
           defaults_t defaultPedestals,
           defaults_t defaultGains,
           defaults_t defaultThresholds_Raw,
           defaults_t defaultThresholds_MeV,
           defaults_t defaultRelativeGains,
           Channel_t::Type_t channelType = Channel_t::Type_t::Integral
           );
    virtual ~Energy();

    const Detector_t::Type_t DetectorType;
    const Channel_t::Type_t ChannelType; // can be Integral or IntegralShort

    const std::shared_ptr<DataManager> calibrationManager;

    const Calibration::Converter::ptr_t Converter;

    CalibType Pedestals;
    CalibType Gains;
    CalibType Thresholds_Raw;
    CalibType Thresholds_MeV;
    CalibType RelativeGains;

    std::vector<CalibType*> AllCalibrations = {
        std::addressof(Pedestals),
        std::addressof(Gains),
        std::addressof(Thresholds_Raw),
        std::addressof(Thresholds_MeV),
        std::addressof(RelativeGains)
    };

};

}}  // namespace ant::calibration
