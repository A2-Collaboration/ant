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

class TAPS_ShortEnergy : public Energy
{


public:
    class ThePhysics : public analysis::Physics {

    protected:
        TH2D* h_pedestals = nullptr;

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
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::Manager_traits> >& guis) override;
protected:
    std::shared_ptr<expconfig::detector::TAPS> taps_detector;
};

}} // namespace ant::calibration
