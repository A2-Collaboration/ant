#pragma once

#include "calibration/Calibration.h"
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
    PID_Energy(
            std::shared_ptr<expconfig::detector::PID> pid,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            const std::vector<double>& defaultPedestals = {100},
            const std::vector<double>& defaultGains = {0.014},
            const std::vector<double>& defaultThresholds = {0.001},
            const std::vector<double>& defaultRelativeGains = {1.0}
    );

    virtual ~PID_Energy();

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;

    using Energy::ApplyTo; // still this ApplyTo should be used
    virtual void ApplyTo(TEventData& reconstructed) override;

protected:
    std::shared_ptr<expconfig::detector::PID> pid_detector;

};

}} // namespace ant::calibration
