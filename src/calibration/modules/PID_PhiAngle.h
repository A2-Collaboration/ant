#pragma once

#include "Calibration.h"

namespace ant {
namespace calibration {

class PID_PhiAngle : public Calibration::PhysicsModule
{
public:
    PID_PhiAngle() :
        PhysicsModule("PID_PhiAngle") {}


    class ThePhysics : public Physics {
    public:
        using Physics::Physics;

        virtual void ProcessEvent(const Event& event) override;
        virtual void Finish() override ;
        virtual void ShowResult() override;
    };

    virtual std::unique_ptr<Physics> GetPhysicsModule() {
        return std_ext::make_unique<ThePhysics>(GetName());
    }
};

}}