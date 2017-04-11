#pragma once

#include "Energy.h"
#include "Energy_GUI.h"

class TH1;
class TH1D;

namespace ant {

class TH2TAPS;

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

    using detector_ptr_t = std::shared_ptr<const expconfig::detector::TAPS>;

    struct GUI_Gains : GUI_CalibType {
        GUI_Gains(const std::string& basename,
                  OptionsPtr options,
                  CalibType& type,
                  const std::shared_ptr<DataManager>& calmgr,
                  const detector_ptr_t& taps);
        virtual ~GUI_Gains();

        virtual void InitGUI(gui::ManagerWindow_traits& window) override;
        virtual DoFitReturn_t DoFit(const TH1& hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;
        // change the histogram name here explicitly to match physics analysis class
        // TAPS_ShortEnergy is just different than CB_Energy or TAPS_Energy
        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
    protected:
        std::shared_ptr<gui::PeakingFitFunction> func;
        gui::CalCanvas* canvas;
        TH1*  h_projection = nullptr;
        TH1D* h_peaks = nullptr;
        TH1D* h_relative = nullptr;

        TH2TAPS* h_peaks_taps = nullptr;
        TH2TAPS* h_relative_taps = nullptr;

        const detector_ptr_t taps_detector;
    };

    struct GUI_Pedestals : energy::GUI_Pedestals {
        GUI_Pedestals(const std::string& basename,
                      OptionsPtr options,
                      CalibType& type,
                      const std::shared_ptr<DataManager>& calmgr,
                      const detector_ptr_t& taps,
                      std::shared_ptr<gui::PeakingFitFunction> fitfunction);
        virtual DoFitReturn_t DoFit(const TH1& hist, unsigned channel) override;
    protected:
        const detector_ptr_t taps_detector;
    };

    TAPS_ShortEnergy(
            const detector_ptr_t& taps,
            const std::shared_ptr<DataManager>& calmgr,
            Calibration::Converter::ptr_t converter,
            /// \todo Do not provide defaults in ctor, setups should decide
            defaults_t defaultPedestals = {100},
            defaults_t defaultGains = {0.3},
            defaults_t defaultThresholds_Raw = {0.0},
            defaults_t defaultThresholds_MeV = {0.0},
            defaults_t defaultRelativeGains = {1.0});

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;
protected:

    detector_ptr_t taps_detector;
};

}} // namespace ant::calibration
