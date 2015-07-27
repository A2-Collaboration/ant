#pragma once

#include "calibration/Calibration.h"
#include "calibration/modules/Integral.h"

class TH1;

namespace ant {
namespace calibration {

class CB_Energy : public Integral
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

    CB_Energy(Calibration::Converter::ptr_t converter,
              const double defaultPedestal,
              const double defaultGain,
              const double defaultThreshold);

    // BaseModule interface
    virtual std::unique_ptr<Physics> GetPhysicsModule();
};

}} // namespace ant::calibration
