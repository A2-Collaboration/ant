#pragma once

#include "Energy.h"

class TH1;

namespace ant {

namespace expconfig {
namespace detector {
struct PID;
}}

namespace calibration {

class PID_Energy :
        public Energy,
        public ReconstructHook::EventData
{


public:

    using detector_ptr_t = std::shared_ptr<const expconfig::detector::PID>;

    PID_Energy(
            const detector_ptr_t& pid,
            const std::shared_ptr<DataManager>& calmgr,
            const Calibration::Converter::ptr_t& converter,
            defaults_t defaultPedestals = {100},
            defaults_t defaultGains = {0.014},
            defaults_t defaultThresholds_Raw = {15},  // roughly width of pedestal peak
            defaults_t defaultThresholds_MeV = {0.1}, // highly depends on calibration!
            defaults_t defaultRelativeGains = {1.0}
    );

    virtual ~PID_Energy();

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;

    using Energy::ApplyTo; // still this ApplyTo should be used
    virtual void ApplyTo(TEventData& reconstructed) override;

protected:
    const detector_ptr_t pid_detector;

};

}} // namespace ant::calibration
