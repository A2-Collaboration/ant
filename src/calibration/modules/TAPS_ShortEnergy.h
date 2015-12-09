#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;
class TH1D;

namespace ant {

namespace expconfig {
namespace detector {
struct TAPS;
}}

namespace calibration {

namespace gui {
class PeakingFitFunction;
}

class TAPS_ShortEnergy : public Energy
{

public:

    /// \todo Add some GUI to calibrate the short gain

    struct GUI_Pedestals : Energy::GUI_Pedestals {
        GUI_Pedestals(const std::string& basename,
                      CalibType& type,
                      const std::shared_ptr<DataManager>& calmgr,
                      const std::shared_ptr<expconfig::detector::TAPS>& taps,
                      std::shared_ptr<gui::PeakingFitFunction> fitfunction);
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) override;
    protected:
        std::shared_ptr<expconfig::detector::TAPS> taps_detector;
    };

    TAPS_ShortEnergy(
            std::shared_ptr<expconfig::detector::TAPS> taps,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            double defaultPedestal = 100,
            double defaultGain = 0.3,
            double defaultThreshold = 0,
            double defaultRelativeGain = 1.0);

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis) override;
protected:


    std::shared_ptr<expconfig::detector::TAPS> taps_detector;
};

}} // namespace ant::calibration
