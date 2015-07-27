#pragma once

#include "calibration/Calibration.h"
#include "Integral.h"

class TH1;

namespace ant {
namespace calibration {

class TAPS_Energy : public Integral
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

    TAPS_Energy(Calibration::Converter::ptr_t converter,
              const double defaultPedestal = 100,
              const double defaultGain = 0.3,
              const double defaultThreshold = 1);

    // BaseModule interface
    virtual std::unique_ptr<Physics> GetPhysicsModule();
};

}} // namespace ant::calibration
