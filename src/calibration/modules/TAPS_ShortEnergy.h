#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {

namespace expconfig {
namespace detector {
class TAPS;
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
    protected:
        std::shared_ptr<gui::PeakingFitFunction> func;
        gui::CalCanvas* canvas;
        TH1*  h_projection = nullptr;
        TH1D* h_peaks = nullptr;
        TH1D* h_relative = nullptr;

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

    class ThePhysics : public analysis::Physics {

    protected:


        TH2D* h_pedestals = nullptr;
        TH2D* h_rel_gamma = nullptr;

        std::shared_ptr<expconfig::detector::TAPS> taps_detector;

    public:



        ThePhysics(const std::string& name,
                   std::shared_ptr<expconfig::detector::TAPS> taps);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override;
        virtual void ShowResult() override;
    };

    TAPS_ShortEnergy(
            std::shared_ptr<expconfig::detector::TAPS> taps,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            double defaultPedestal = 100,
            double defaultGain = 0.3,
            double defaultThreshold = 0,
            double defaultRelativeGain = 1.0);

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis) override;
protected:


    std::shared_ptr<expconfig::detector::TAPS> taps_detector;
};

}} // namespace ant::calibration
