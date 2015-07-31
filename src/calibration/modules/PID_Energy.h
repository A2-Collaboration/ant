#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {
namespace calibration {

class PID_Energy : public Energy
{


public:
    class ThePhysics : public Physics {

    protected:
        TH2* ggIM = nullptr;

    public:
        ThePhysics(const std::string& name);

        virtual void ProcessEvent(const Event& event);
        virtual void Finish();
        virtual void ShowResult();
    };

    PID_Energy(std::shared_ptr<CalibrationDataManager> calmgr,
               Calibration::Converter::ptr_t converter,
               double defaultPedestal = 100,
               double defaultGain = 0.014,
               double defaultThreshold = 0.001,
               double defaultRelativeGain = 1.0);

    // BaseModule interface
    virtual std::unique_ptr<Physics> GetPhysicsModule();
};

}} // namespace ant::calibration
