#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {

namespace expconfig {
namespace detector {
class CB;
}}

namespace calibration {

class CB_Energy : public Energy
{


public:
    struct TheGUI : GUI_CalibType {
        TheGUI(const std::string& basename, CalibType& type, CB_Energy* parent);

        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI() override;
        virtual std::list<gui::CalCanvas*> GetCanvases() const override;
        virtual bool DoFit(TH1* hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishRange() override;
    protected:
        CB_Energy* p;
        gui::CalCanvas* c_fit;
        gui::CalCanvas* c_overview;
    };
    friend class TheGUI;

    struct ThePhysics : Physics {

        ThePhysics(const std::string& name, const std::string& hist_name, unsigned nChannels);

        virtual void ProcessEvent(const Event& event) override;
        virtual void Finish() override;
        virtual void ShowResult() override;

    protected:
        TH2* ggIM = nullptr;
    };

    CB_Energy(
            std::shared_ptr<expconfig::detector::CB> cb,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            double defaultPedestal = 0,
            double defaultGain = 0.07,
            double defaultThreshold = 2,
            double defaultRelativeGain = 1.0);

    virtual std::unique_ptr<Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis) override;


protected:
    std::shared_ptr<expconfig::detector::CB> cb_detector;
};

}} // namespace ant::calibration
