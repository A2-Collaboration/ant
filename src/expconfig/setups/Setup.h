#pragma once

#include "ExpConfig.h"
#include "SetupRegistry.h"

#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"

#include "tree/THeaderInfo.h"

#include "base/Paths.h"
#include "base/std_ext.h"
#include "base/Logger.h"

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
#include "calibration/modules/CB_TimeWalk.h"
#include "calibration/modules/PID_Energy.h"
#include "calibration/modules/PID_PhiAngle.h"
#include "calibration/modules/TAPS_Energy.h"
#include "calibration/modules/TAPS_ShowerCorrection.h"
#include "calibration/modules/TAPSVeto_Energy.h"

#include "calibration/converters/CATCH_TDC.h"
#include "calibration/converters/MultiHit16bit.h"
#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/ScalerFrequency.h"

#include "calibration/DataManager.h"

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

    virtual std::list< std::shared_ptr<Updateable_traits> > GetUpdateables() const override {
        std::list< std::shared_ptr<Updateable_traits> > list;
        for(const auto& hook : reconstruct_hooks) {
            std_ext::AddToSharedPtrList<Updateable_traits, ReconstructHook::Base>(
                        hook, list
                        );
        }
        for(const auto& calib : calibrations) {
            std_ext::AddToSharedPtrList<Updateable_traits, Calibration::BaseModule>(
                        calib, list
                        );
        }
        for(const auto& det : detectors) {
            std_ext::AddToSharedPtrList<Updateable_traits, Detector_t>(
                        det, list
                        );
        }
        return list;
    }


    virtual std::string GetPIDCutsDirectory() const override {
        return std::string(ANT_PATH_DATABASE)+"/"+GetName()+"/cuts";
    }

    std::string GetName() const override final {
        return name_;
    }

protected:
    Setup(const std::string& name) :
        name_(name)
    {
        std::string filename = std::string(ANT_PATH_DATABASE)+"/"+GetName()+"/calibration.root";
        calibrationDataManager = std::make_shared<calibration::DataManager>(filename);
    }

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
        for(auto detector : detectors) {
            auto cfg = std::dynamic_pointer_cast<UnpackerAcquConfig, Detector_t>(detector);
            if(cfg == nullptr)
                continue;
            //std::vector<hit_mapping_t> hit_mappings_;
            //std::vector<scaler_mapping_t> scaler_mappings_;
            cfg->BuildMappings(hit_mappings, scaler_mappings);
            /// \todo check that the detectors do not add overlapping mappings
        }
    }

    void IgnoreDetectorChannel(Detector_t::Type_t type, unsigned channel) {
        for(auto& detector : detectors) {
            if(detector->Type == type) {
                detector->SetIgnored(channel);
                return;
            }
        }
        LOG(WARNING) << "Setup " << GetName() << " ignored non-existing channel "
                     << channel << " of detector " << Detector_t::ToString(type);
    }
    void IgnoreDetectorChannels(Detector_t::Type_t type, const std::vector<unsigned>& channels) {
        for(unsigned channel : channels)
            IgnoreDetectorChannel(type, channel);
    }

    static std::shared_ptr<calibration::DataManager> CreateCalibrationDataManager(const std::string& setupname) {
        std::string filename = std::string(ANT_PATH_DATABASE)+"/"+setupname+"/calibration.root";
        return std::make_shared<calibration::DataManager>(filename);
    }

    std::list< std::shared_ptr<Detector_t> > detectors;
    std::list< std::shared_ptr<ReconstructHook::Base> > reconstruct_hooks;
    std::list< std::shared_ptr<Calibration::BaseModule> > calibrations;

    std::shared_ptr<calibration::DataManager> calibrationDataManager;
private:
    const std::string name_;
};

}} // namespace ant::expconfig
