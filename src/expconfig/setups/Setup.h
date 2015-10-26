#pragma once

#include "ExpConfig.h"
#include "SetupRegistry.h"

#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"

#include "tree/THeaderInfo.h"

#include "calibration/Calibration.h"
#include "calibration/DataManager.h"


#include <functional>
#include <string>

namespace ant {
namespace expconfig {

using SetupOptPtr = std::shared_ptr<const OptionsList>;


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
private:
    const std::string name_;
    SetupOptPtr options;

public:
    virtual std::list< std::shared_ptr< Calibration::PhysicsModule> > GetCalibrations() const override;

    virtual std::list< std::shared_ptr< ReconstructHook::Base > > GetReconstructHooks() const override;

    virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const override {
        return detectors;
    }

    virtual std::list< std::shared_ptr<Updateable_traits> > GetUpdateables() const override;


    virtual std::string GetPIDCutsDirectory() const override;

    virtual std::shared_ptr<calibration::DataManager> GetCalibrationDataManager() const override final {
        return calibrationDataManager;
    }

    virtual std::string GetName() const override final {
        return name_;
    }

protected:
    Setup(const std::string& name, SetupOptPtr opt);

    std::string GetOption(const std::string& key) const;
    bool IsFlagSet(const std::string& key) const;

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

    bool Matches(const THeaderInfo& header) const override;

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const;

    void IgnoreDetectorChannel(Detector_t::Type_t type, unsigned channel);
    void IgnoreDetectorChannels(Detector_t::Type_t type, const std::vector<unsigned>& channels);

    std::list< std::shared_ptr<Detector_t> > detectors;
    std::list< std::shared_ptr<ReconstructHook::Base> > reconstruct_hooks;
    std::list< std::shared_ptr<Calibration::BaseModule> > calibrations;

    std::shared_ptr<calibration::DataManager> calibrationDataManager;

};

}} // namespace ant::expconfig
