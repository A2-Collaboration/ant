#pragma once

#include "expconfig/ExpConfig.h"

#include "SetupRegistry.h"

#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"

#include "tree/TID.h"

#include "calibration/Calibration.h"
#include "calibration/DataManager.h"

#include "base/OptionsList.h"

#include <functional>
#include <string>

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
        public UnpackerAcquConfig,
        public UnpackerA2GeantConfig
{
private:
    const std::string name_;

public:
    virtual std::list< std::shared_ptr< Calibration::PhysicsModule> > GetCalibrations() const override;

    virtual std::list< std::shared_ptr< ReconstructHook::Base > > GetReconstructHooks() const override;

    virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const override {
        return detectors;
    }

    virtual std::list< std::shared_ptr<Updateable_traits> > GetUpdateables() const override;


    virtual std::string GetPIDCutsDirectory() const override;
    virtual std::string GetPhysicsFilesDirectory() const override;

    virtual std::shared_ptr<calibration::DataManager> GetCalibrationDataManager() const override final {
        return calibrationDataManager;
    }

    virtual std::string GetName() const override final {
        return name_;
    }

    virtual double GetElectronBeamEnergy() const override {
        return std::numeric_limits<double>::quiet_NaN();
    }

protected:
    const OptionsPtr Options;

    Setup(const std::string& name, OptionsPtr opt);

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

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const override;

    void IgnoreDetectorChannel(Detector_t::Type_t type, unsigned channel);
    void IgnoreDetectorChannels(Detector_t::Type_t type, const std::vector<unsigned>& channels);

    std::list< std::shared_ptr<Detector_t> > detectors;
    std::list< std::shared_ptr<ReconstructHook::Base> > reconstruct_hooks;
    std::list< std::shared_ptr<Calibration::BaseModule> > calibrations;

    std::shared_ptr<calibration::DataManager> calibrationDataManager;

};

}} // namespace ant::expconfig
