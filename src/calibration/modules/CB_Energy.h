#pragma once

#include "calibration/Calibration.h"
#include "Energy.h"

class TH1;

namespace ant {
namespace calibration {

class CB_Energy : public Energy
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
              const double defaultPedestal = 0,
              const double defaultGain = 0.07,
              const double defaultThreshold = 2);

    // BaseModule interface
    virtual std::unique_ptr<Physics> GetPhysicsModule();
};

}} // namespace ant::calibration
