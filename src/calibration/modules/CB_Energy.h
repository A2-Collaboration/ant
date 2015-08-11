#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {

class TH2CB;

namespace expconfig {
namespace detector {
class CB;
}}

namespace calibration {

namespace gui {
class FitGausPol3;
}

class CB_Energy : public Energy
{


public:
    struct TheGUI : GUI_CalibType {
        TheGUI(const std::string& basename,
               CalibType& type,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<expconfig::detector::CB>& cb);

        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI() override;
        virtual std::list<gui::CalCanvas*> GetCanvases() const override;
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishRange() override;
        virtual void StoreFinishRange(const interval<TID>& range) override;
    protected:
        std::shared_ptr<expconfig::detector::CB> cb_detector;
        std::shared_ptr<gui::FitGausPol3> func;
        TH1* projection = nullptr;
        gui::CalCanvas* c_fit;
        gui::CalCanvas* c_overview;
        TH2CB* hist_cb;
    };

    struct ThePhysics : analysis::Physics {

        ThePhysics(const std::string& name, const std::string& hist_name, unsigned nChannels);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
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

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis) override;


protected:
    std::shared_ptr<expconfig::detector::CB> cb_detector;
};

}} // namespace ant::calibration
