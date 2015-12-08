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
    };

//    struct ThePhysics : analysis::Physics {

//        ThePhysics(const std::string& name, const std::string& hist_name, unsigned nChannels);

//        virtual void ProcessEvent(const analysis::data::Event& event) override;
//        virtual void Finish() override;
//        virtual void ShowResult() override;

//    protected:
//        TH2* ggIM = nullptr;
//    };

    CB_Energy(
            std::shared_ptr<expconfig::detector::CB> cb,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            double defaultPedestal = 0,
            double defaultGain = 0.07,
            double defaultThreshold = 2,
            double defaultRelativeGain = 1.0);

    virtual std::vector<std::string> GetPhysicsModules() const override;
    virtual void GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis) override;


protected:
    std::shared_ptr<expconfig::detector::CB> cb_detector;

    // CB has online pedestal subtraction
    virtual bool NeedsPedestals() const override { return false; }
};

}} // namespace ant::calibration
