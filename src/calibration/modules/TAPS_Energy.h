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

class TAPS_Energy : public Energy
{


public:
    class ThePhysics : public Physics {

    protected:
        TH2* ggIM = nullptr;
        std::shared_ptr<expconfig::detector::TAPS> taps_detector;

    public:
        ThePhysics(const std::string& name,
                   std::shared_ptr<expconfig::detector::TAPS> taps);

        virtual void ProcessEvent(const Event& event);
        virtual void Finish();
        virtual void ShowResult();
    };

    TAPS_Energy(
            std::shared_ptr<expconfig::detector::TAPS> taps,
            std::shared_ptr<DataManager> calmgr,
            Calibration::Converter::ptr_t converter,
            double defaultPedestal = 100,
            double defaultGain = 0.3,
            double defaultThreshold = 1,
            double defaultRelativeGain = 1.0);

    virtual std::unique_ptr<Physics> GetPhysicsModule();
    virtual std::list<std::unique_ptr<calibration::gui::Manager_traits>> GetGUIs() override {
        return {};
    }
protected:
    std::shared_ptr<expconfig::detector::TAPS> taps_detector;
};

}} // namespace ant::calibration
