#pragma once

#include "ExpConfig.h"
#include "unpacker/UnpackerAcqu.h"
#include "base/std_ext.h"
#include "tree/THeaderInfo.h"

// this header always includes all detectors and calibrations
// for convinient access in the derived classes

#include "detectors/Trigger.h"
#include "detectors/CB.h"
#include "detectors/PID.h"
#include "detectors/TAPS.h"
#include "detectors/TAPSVeto.h"
#include "detectors/EPT.h"


#include "calibration/modules/EnergyInvariantMass.h"
#include "calibration/modules/IntegralSADC.h"
#include "calibration/modules/IntegralTAPS.h"
#include "calibration/modules/IntegralCAEN.h"
#include "calibration/modules/TimingCATCH.h"
#include "calibration/modules/TimingTAPS.h"

#include <functional>

namespace ant {
namespace expconfig {

class Setup :
        public ExpConfig::Module,
        public ExpConfig::Reconstruct,
        public UnpackerAcquConfig
{
public:
    // every setup should have a name,
    // does not need to be unique as long as the matching
    // is unique
    /// \todo Implement cmdline match override to make this meaningful...
    virtual std::string GetName() const = 0;

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

// stuff for semiauto-registering the setups

template<class T>
std::shared_ptr<Setup> setup_factory()
{
    return std::move(std::make_shared<T>());
}

using setup_creator = std::function<std::shared_ptr<Setup>()>;

class SetupRegistry
{
private:
    using setup_creators_t = std::list<setup_creator>;
    using setups_t = std::list< std::shared_ptr<Setup> >;
    setup_creators_t setup_creators;
    setups_t setups;
    void init_setups();
public:
    static SetupRegistry& get();
    void add(setup_creator);
    setups_t::iterator begin();
    setups_t::iterator end();
};

class SetupRegistration
{
public:
    SetupRegistration(setup_creator);
};

#define AUTO_REGISTER_SETUP(setup) \
    SetupRegistration _setup_registration_ ## setup(setup_factory<setup>);

}} // namespace ant::expconfig
