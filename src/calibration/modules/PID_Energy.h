#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {

namespace expconfig {
namespace detector {
class PID;
}}

namespace calibration {

namespace gui {
class FitGausPol0;
}

class PID_Energy : public Energy
{


public:
    PID_Energy(
            std::shared_ptr<expconfig::detector::PID> pid,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            double defaultPedestal = 100,
            double defaultGain = 0.014,
            double defaultThreshold = 0.001,
            double defaultRelativeGain = 1.0);

    virtual ~PID_Energy();


    class ThePhysics : public analysis::Physics {

    protected:
        TH2D* h_pedestals = nullptr;
        TH3D* h_bananas = nullptr;

    public:
         ThePhysics(const std::string& name, unsigned nChannels);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override;
        virtual void ShowResult() override;
    }; // ThePhysics

    struct GUI_Pedestals : GUI_CalibType {
        GUI_Pedestals(const std::string& basename,
               CalibType& type,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<Detector_t>& detector);

        virtual void InitGUI(gui::ManagerWindow_traits* window) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishRange() override;
    protected:
        std::shared_ptr<gui::FitGausPol0> func;
        gui::CalCanvas* canvas;
        TH1*  h_projection = nullptr;

    }; // TheGUI

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::Manager_traits> >& guis) override;

protected:
    std::shared_ptr<expconfig::detector::PID> pid_detector;

};

}} // namespace ant::calibration
