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
    struct ThePhysics : Physics {
        ThePhysics(const std::string& name, unsigned nChannels);

        virtual void ProcessEvent(const Event& event);
        virtual void Finish();
        virtual void ShowResult();

    protected:
        TH2* ggIM = nullptr;

    };

    struct TheGUI : gui::Manager_traits {

        virtual std::string GetHistogramName() const;
        virtual unsigned GetNumberOfChannels() const;
        virtual void InitGUI();
        virtual std::list<gui::CalCanvas*> GetCanvases() const;
        virtual void StartRange(const interval<TID>& interval);
        virtual bool DoFit(TH1* hist, unsigned channel);
        virtual void DisplayFit();
        virtual void StoreFit(unsigned channel);
        virtual bool FinishRange();
        virtual void StoreFinishRange(const interval<TID>& interval);

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
    virtual std::list<std::unique_ptr<calibration::gui::Manager_traits>> GetGUIs() override;


protected:
    std::shared_ptr<expconfig::detector::CB> cb_detector;
};

}} // namespace ant::calibration
