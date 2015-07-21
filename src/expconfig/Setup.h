#ifndef ANT_EXPCONFIG_SETUP
#define ANT_EXPCONFIG_SETUP

#include "ExpConfig.h"
#include "unpacker/UnpackerAcqu.h"
#include "base/std_ext.h"

// this header always includes all detectors and calibrations
// for convinient access in the derived classes

#include "detectors/Trigger.h"
#include "detectors/CB.h"
#include "detectors/PID.h"
#include "detectors/TAPS.h"
#include "detectors/EPT.h"

#include "calibration/modules/EnergyInvariantMass.h"
#include "calibration/modules/IntegralSADC.h"
#include "calibration/modules/IntegralTAPS.h"
#include "calibration/modules/IntegralCAEN.h"
#include "calibration/modules/TimingCATCH.h"
#include "calibration/modules/TimingTAPS.h"


namespace ant {
namespace expconfig {

class Setup :
        public ExpConfig::Module,
        public ExpConfig::Reconstruct,
        public UnpackerAcquConfig
{

protected:
    Setup() = default;

    void AddDetector(const std::shared_ptr<Detector_t>& detector) {
        detectors.push_back(detector);
    }

    template<typename T, typename... Args>
    void AddDetector(Args&&... args) {
        AddDetector(std::make_shared<T>(std::forward<Args>(args)...));
    }

    template<typename T, typename... Args>
    void AddCalibration(Args&&... args) {
        calibrations.push_back(std::make_shared<T>(std::forward<Args>(args)...));
    }

    virtual std::list< std::shared_ptr< CalibrationApply_traits > > GetCalibrations() const override {
        return calibrations;
    }

    virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const override {
        return detectors;
    }

    bool Matches(const THeaderInfo& header) const override {
        // check that all detectors match
        for(const auto& detector : detectors) {
            const auto& ptr = dynamic_pointer_cast<ExpConfig::Base, Detector_t>(detector);
            if(ptr == nullptr)
                continue;
            if(!ptr->Matches(header))
                return false;
        }
        return true;
    }

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const
    {
        // the base setup simply asks its underlying
        // detectors for the mappings
        for(const auto& detector : detectors) {
            const UnpackerAcquConfig* cfg
                    = dynamic_cast<const UnpackerAcquConfig*>(detector.get());
            if(cfg == nullptr)
                continue;
            cfg->BuildMappings(hit_mappings, scaler_mappings);
        }
    }

    std::list< std::shared_ptr<Detector_t> > detectors;
    std::list< std::shared_ptr<CalibrationApply_traits> > calibrations;
};

}} // namespace ant::expconfig

#endif
