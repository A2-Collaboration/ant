#pragma once

#include "ExpConfig.h"

#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"

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

#include "calibration/modules/Time.h"
#include "calibration/modules/Scaler.h"
#include "calibration/modules/TAPS_ShowerCorrection.h"
#include "calibration/modules/CB_Energy.h"
#include "calibration/modules/PID_Energy.h"
#include "calibration/modules/TAPS_Energy.h"
#include "calibration/modules/TAPSVeto_Energy.h"
#include "calibration/CalibrationDataManager.h"


#include "calibration/converters/CATCH_TDC.h"
#include "calibration/converters/MultiHit16bit.h"
#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/ScalerFrequency.h"



#include <functional>

namespace ant {
namespace expconfig {

class Setup :
        public ExpConfig::Setup,
        public ExpConfig::Reconstruct,
        public UnpackerAcquConfig,
        public UnpackerA2GeantConfig
{
public:
    virtual std::list< std::shared_ptr< Calibration::PhysicsModule> > GetCalibrations() override {
        // search the hooks for modules which are physics modules
        std::list< std::shared_ptr<Calibration::PhysicsModule> > calibrations;
        for(const auto& hook : reconstruct_hooks) {
            std_ext::AddToSharedPtrList<Calibration::PhysicsModule, ReconstructHook::Base>(
                        hook, calibrations
                        );
        }
       return calibrations;
    }

protected:
    Setup() = default;

    void AddDetector(const std::shared_ptr<Detector_t>& detector) {
        detectors.push_back(detector);
    }
    template<typename T, typename... Args>
    void AddDetector(Args&&... args) {
        AddDetector(std::make_shared<T>(std::forward<Args>(args)...));
    }

    void AddHook(const std::shared_ptr<ReconstructHook::Base>& hook) {
        reconstruct_hooks.push_back(hook);
    }
    template<typename T, typename... Args>
    void AddHook(Args&&... args) {
        AddHook(std::make_shared<T>(std::forward<Args>(args)...));
    }

    virtual std::list< std::shared_ptr< ReconstructHook::Base > > GetReconstructHooks() const override {
        return reconstruct_hooks;
    }

    virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const override {
        return detectors;
    }

    bool Matches(const THeaderInfo& header) const override {
        // check that all detectors match
        for(const auto& detector : detectors) {
            const auto& ptr = std::dynamic_pointer_cast<ExpConfig::Base, Detector_t>(detector);
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
            //std::vector<hit_mapping_t> hit_mappings_;
            //std::vector<scaler_mapping_t> scaler_mappings_;
            cfg->BuildMappings(hit_mappings, scaler_mappings);
            /// \todo check that the detectors do not add overlapping mappings
        }
    }

    std::list< std::shared_ptr<Detector_t> > detectors;
    std::list< std::shared_ptr<ReconstructHook::Base> > reconstruct_hooks;
};

// stuff for semiauto-registering the setups
// don't forget to use AUTO_REGISTER_SETUP

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
