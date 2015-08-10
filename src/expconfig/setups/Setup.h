#pragma once

#include "ExpConfig.h"
#include "SetupRegistry.h"

#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"

#include "base/Paths.h"
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
#include "calibration/modules/CB_Energy.h"
#include "calibration/modules/PID_Energy.h"
#include "calibration/modules/PID_PhiAngle.h"
#include "calibration/modules/TAPS_Energy.h"
#include "calibration/modules/TAPS_ShowerCorrection.h"
#include "calibration/modules/TAPSVeto_Energy.h"

#include "calibration/converters/CATCH_TDC.h"
#include "calibration/converters/MultiHit16bit.h"
#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/ScalerFrequency.h"

#include "calibration/CalibrationDataManager.h"

#include <functional>

namespace ant {
namespace expconfig {


/**
 * @brief The Setup class serves as a base class for all known beamtime setup configurations
 *
 * The class provides useful routines and stores the instantiated calibrations and detectors,
 * which are used to automatically construct the mappings for the unpacker, for example.
 *
 */
class Setup :
        public ExpConfig::Setup,
        public ExpConfig::Reconstruct,
        public UnpackerAcquConfig,
        public UnpackerA2GeantConfig
{
public:
    virtual std::list< std::shared_ptr< Calibration::PhysicsModule> > GetCalibrations() const override {
        // search the hooks for modules which are physics modules
        std::list< std::shared_ptr<Calibration::PhysicsModule> > list;
        for(const auto& hook : reconstruct_hooks) {
            std_ext::AddToSharedPtrList<Calibration::PhysicsModule, ReconstructHook::Base>(
                        hook, list
                        );
        }
        for(const auto& calib : calibrations) {
            std_ext::AddToSharedPtrList<Calibration::PhysicsModule, Calibration::BaseModule>(
                        calib, list
                        );
        }
        return list;
    }

    virtual std::list< std::shared_ptr< ReconstructHook::Base > > GetReconstructHooks() const override {
        std::list< std::shared_ptr< ReconstructHook::Base > > list = reconstruct_hooks;
        for(const auto& calib : calibrations) {
            std_ext::AddToSharedPtrList<ReconstructHook::Base, Calibration::BaseModule>(
                        calib, list
                        );
        }
        return list;
    }

    virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const override {
        return detectors;
    }

    virtual std::string GetPIDCutsDirectory() const override {
        return std::string(ANT_PATH_DATABASE)+"/"+GetName()+"/cuts";
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

    void AddCalibration(const std::shared_ptr<Calibration::BaseModule>& calib) {
        calibrations.push_back(calib);
    }
    template<typename T, typename... Args>
    void AddCalibration(Args&&... args) {
        AddCalibration(std::make_shared<T>(std::forward<Args>(args)...));
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
    std::list< std::shared_ptr<Calibration::BaseModule> > calibrations;
};

}} // namespace ant::expconfig
