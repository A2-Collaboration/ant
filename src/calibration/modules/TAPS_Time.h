#pragma once

#include "Time.h"


namespace ant {

namespace expconfig {
namespace detector {
class TAPS;
}}

namespace calibration {

class TAPS_Time : public Time
{

public:
    TAPS_Time(std::shared_ptr<expconfig::detector::TAPS> taps,
              std::shared_ptr<DataManager> calmgr,
              Calibration::Converter::ptr_t converter,
              const interval<double>& timeWindow_BaF2 = {-std_ext::inf, std_ext::inf},
              const interval<double>& timeWindow_PbWO4 = {-std_ext::inf, std_ext::inf});

    class ThePhysics : public Time::ThePhysics {
    public:
        ThePhysics(const std::string& name, const std::string& histName,
                   const std::shared_ptr<Detector_t>& theDetector);
        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void ShowResult() override;
    protected:
        TH2D* hTimeToTagger;
    }; // ThePhysics

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;


};

}}  // namespace ant::calibration
