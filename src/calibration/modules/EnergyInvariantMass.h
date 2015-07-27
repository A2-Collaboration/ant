#pragma once

#include "Calibration.h"

class TH1;

namespace ant {
namespace calibration {

class EnergyInvariantMass :
        public Calibration::PhysicsModule,
        public ReconstructHook::DetectorReadHits
{


public:
    class ThePhysics : public Physics {

    protected:
        TH1* ggIM = nullptr;

    public:
        ThePhysics(const std::string& name);

        virtual void ProcessEvent(const Event& event);
        virtual void Finish();
        virtual void ShowResult();
    };

    EnergyInvariantMass();

    // ReconstructHook
    virtual void ApplyTo(const readhits_t&, extrahits_t&) override {}

    // CalibrationUpdate_traits interface
    virtual std::list<TID> GetChangePoints() const override { return {}; }
    virtual void Update(const TID &) override {}

    // PhysicsFactory interface
    virtual std::unique_ptr<Physics> GetPhysicsModule();
};

}} // namespace ant::calibration
