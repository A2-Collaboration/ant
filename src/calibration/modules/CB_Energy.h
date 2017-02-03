#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1D;

namespace ant {

class TH2CB;

namespace expconfig {
namespace detector {
struct CB;
}}

namespace calibration {

namespace gui {
class FitGausPol3;
}

class CB_Energy : public Energy
{


public:
    struct GUI_Gains : GUI_CalibType {
        GUI_Gains(const std::string& basename,
               OptionsPtr options,
               CalibType& type,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<Detector_t>& detector);

        virtual void InitGUI(gui::ManagerWindow_traits* window) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;
    protected:
        std::shared_ptr<gui::FitGausPol3> func;
        gui::CalCanvas* canvas;
        TH1*  h_projection = nullptr;
        TH1D* h_peaks = nullptr;
        TH1D* h_relative = nullptr;

        TH2CB* h_peaks_cb = nullptr;
        TH2CB* h_relative_cb = nullptr;

        double AutoStopOnChi2 = 6;
        interval<double> FitRange = {20, 200};
        double ConvergenceFactor = 1.0;
    };

    CB_Energy(
            std::shared_ptr<expconfig::detector::CB> cb,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            const std::vector<double>& defaultPedestals,
            const std::vector<double>& defaultGains,
            const std::vector<double>& defaultThresholds,
            const std::vector<double>& defaultRelativeGains);

    virtual void GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;


protected:
    std::shared_ptr<expconfig::detector::CB> cb_detector;
};

}} // namespace ant::calibration
