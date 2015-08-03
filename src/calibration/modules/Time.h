#pragma once

#include "calibration/Calibration.h"
#include "expconfig/Detector_t.h"

#include "tree/TDataRecord.h" // for TKeyValue, TID
#include "base/interval.h"

#include <memory>
#include <limits>

class TH1;

namespace ant {
namespace calibration {

class Time :
        public Calibration::Module,
        public ReconstructHook::DetectorReadHits
{

public:

    Time(
            Detector_t::Type_t DetectorType,
            Calibration::Converter::ptr_t converter,
            double defaultOffset,
            const interval<double>& timeWindow = {-std_ext::inf, std_ext::inf},
            double defaultGain = 1.0, // default gain is 1.0
            const std::vector< TKeyValue<double> >& gains = {}
            );

    // ReconstructHook
    virtual void ApplyTo(const readhits_t& hits, extrahits_t&) override;

    // Updateable_traits interface
    virtual std::vector<std::list<TID>> GetChangePoints() const override;
    void Update(std::size_t index, const TID&) override;

    class ThePhysics : public Physics {
    public:
        using Physics::Physics;

        virtual void ProcessEvent(const Event& event) override;
        virtual void Finish() override;
        virtual void ShowResult() override;
    };

    // Physics_traits interface
    virtual std::unique_ptr<Physics> GetPhysicsModule() {
        return std_ext::make_unique<ThePhysics>(GetName());
    }

protected:

    const Detector_t::Type_t DetectorType;
    const Calibration::Converter::ptr_t Converter;

    const interval<double> TimeWindow;

    const double DefaultOffset;
    std::vector<double> Offsets;

    const double DefaultGain;
    std::vector<double> Gains;

};

}}  // namespace ant::calibration
