#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

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

    struct GUI_Gains : GUI_CalibType {
        GUI_Gains(const std::string& basename,
                  CalibType& type,
                  const std::shared_ptr<DataManager>& calmgr,
                  const std::shared_ptr<expconfig::detector::TAPS>& taps);
        virtual ~GUI_Gains();

        virtual void InitGUI(gui::ManagerWindow_traits* window) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) override;
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

        std::shared_ptr<expconfig::detector::TAPS> taps_detector;
    };

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
            const std::vector<double>& defaultPedestals = {100},
            const std::vector<double>& defaultGains = {0.3},
            const std::vector<double>& defaultThresholds = {0.0},
            const std::vector<double>& defaultRelativeGains = {1.0});

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;
protected:


    std::shared_ptr<expconfig::detector::TAPS> taps_detector;
};

}} // namespace ant::calibration
