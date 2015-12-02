#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {

namespace expconfig {
namespace detector {
struct TAPS;
}}

namespace calibration {

namespace gui {
class FitGausPol3;
}


class TAPS_Energy : public Energy
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
    };

    class ThePhysics : public analysis::Physics {

    protected:
        TH2D* ggIM = nullptr;
        TH2D* timing_cuts = nullptr;
        TH2D* h_pedestals = nullptr;

        std::shared_ptr<expconfig::detector::TAPS> taps_detector;

    public:

        ThePhysics(const std::string& name,
                   std::shared_ptr<expconfig::detector::TAPS> taps);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override;
        virtual void ShowResult() override;
    };

    TAPS_Energy(
            std::shared_ptr<expconfig::detector::TAPS> taps,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            double defaultPedestal = 100,
            double defaultGain = 0.3,
            double defaultThreshold = 1,
            double defaultRelativeGain = 1.0);

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis) override;
protected:
    std::shared_ptr<expconfig::detector::TAPS> taps_detector;
};

}} // namespace ant::calibration
