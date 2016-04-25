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
class FitGausexpo;
}

class CB_SourceCalib : public Energy
{

public:
    struct GUI_Gains : GUI_CalibType {
        GUI_Gains(const std::string& basename,
               CalibType& type,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<Detector_t>& detector);

        virtual void InitGUI(gui::ManagerWindow_traits* window) override;
        virtual DoFitReturn_t DoFit(TH1 *hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual std::shared_ptr<TH1>  GetHistogram(const WrapTFile &file) const;
        virtual bool FinishSlice() override;
        virtual std::string GetName() const override;
    protected:
        std::shared_ptr<gui::FitGausexpo> func;
        gui::CalCanvas* canvas;
        TH1* h_projection = nullptr;
        TH1D* h_peaks = nullptr;


        double AutoStopOnChi2 = 6;

    };

    CB_SourceCalib(
            std::shared_ptr<expconfig::detector::CB> cb,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            double defaultPedestal,
            double defaultGain,
            double defaultThreshold,
            double defaultRelativeGain);

    virtual void GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis) override;

protected:
     std::shared_ptr<expconfig::detector::CB> cb_detector;
};

}}
