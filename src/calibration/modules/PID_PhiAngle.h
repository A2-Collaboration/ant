#pragma once

#include "Calibration.h"
#include "base/interval.h"

class TGraph;

namespace ant {

namespace expconfig {
namespace detector {
class PID;
}}

namespace calibration {

class DataManager;

namespace gui {
class FitGaus;
}


class PID_PhiAngle : public Calibration::Module
{
public:
    PID_PhiAngle(const std::shared_ptr<expconfig::detector::PID>&  pid,
            const std::shared_ptr<DataManager>& calmgr
            );
    virtual ~PID_PhiAngle();

    class ThePhysics : public analysis::Physics {
    protected:
        TH2* pid_cb_phi_corr;
        const ant::interval<double> theta_range;
    public:
        ThePhysics(const std::string& name, unsigned nChannels);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override ;
        virtual void ShowResult() override;
    }; // ThePhysics

    class TheGUI : public gui::CalibModule_traits {
    protected:
        std::shared_ptr<DataManager> calibrationManager;
        std::shared_ptr<expconfig::detector::PID> pid_detector;
        class _FitGauss;
        std::shared_ptr<_FitGauss> func;

        gui::CalCanvas* canvas;

        TH1*  h_projection = nullptr;
        TGraph* h_result;

        std::vector<double> angles;
        std::vector<double> previousAngles;
        std::map< unsigned, std::vector<double> > fitParameters;

        double phi_offset = std::numeric_limits<double>::quiet_NaN();

        bool IgnorePreviousFitParameters = false;
    public:
        TheGUI(const std::string& basename,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<expconfig::detector::PID>& pid
               );
        virtual ~TheGUI();

        virtual std::string GetHistogramName() const override;
        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI(gui::ManagerWindow_traits*) override;

        virtual void StartSlice(const interval<TID>& range) override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;
        virtual void StoreFinishSlice(const interval<TID>& range) override;
    }; // TheGUI

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis) override;

    // Updateable_traits interface
    virtual std::list<UpdateableItemPtr> GetItems() const override;

protected:
    std::shared_ptr<expconfig::detector::PID> pid_detector;
    std::shared_ptr<DataManager> calibrationManager;



};

}}
