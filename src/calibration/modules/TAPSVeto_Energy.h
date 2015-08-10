#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {
namespace calibration {

class TAPSVeto_Energy : public Energy
{


public:
    class ThePhysics : public analysis::Physics {

    protected:
        TH2* ggIM = nullptr;

    public:
        ThePhysics(const std::string& name);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override;
        virtual void ShowResult() override;
    };

    TAPSVeto_Energy(std::shared_ptr<DataManager> calmgr,
                    Calibration::Converter::ptr_t converter,
                    double defaultPedestal = 100,
                    double defaultGain = 0.010,
                    double defaultThreshold = 0.1,
                    double defaultRelativeGain = 1.0);

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::Manager_traits> >& guis) override {
        guis.clear();
    }
};

}} // namespace ant::calibration
