#pragma once

#include "expconfig/ExpConfig.h"

#include "SetupRegistry.h"

#include "unpacker/UnpackerAcqu.h"
#include "unpacker/UnpackerA2Geant.h"

#include "tree/TID.h"

#include "calibration/Calibration.h"
#include "calibration/DataManager.h"
#include "base/piecewise_interval.h"
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
        public expconfig::Setup_traits,
        public UnpackerAcquConfig,
        public UnpackerA2GeantConfig
{
private:
    const std::string name_;
    const bool includeIgnoredElements;

    ant::PiecewiseInterval<double> prompt = {};
    ant::PiecewiseInterval<double> random = {};

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

    virtual bool GetIncludeIgnoredElements() const override {
        return includeIgnoredElements;
    }

    virtual ant::PiecewiseInterval<double> GetPromptWindows() const override {
        return prompt;
    }
    virtual ant::PiecewiseInterval<double> GetRandomWindows() const override {
        return random;
    }

protected:

    Setup(const std::string& name, OptionsPtr opts);

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

    std::list< std::shared_ptr<Detector_t> > detectors;
    std::list< std::shared_ptr<ReconstructHook::Base> > reconstruct_hooks;
    std::list< std::shared_ptr<Calibration::BaseModule> > calibrations;

    std::shared_ptr<calibration::DataManager> calibrationDataManager;

    void AddPromptRange(const interval<double>& i) {
        prompt.emplace_back(i);
    }

    void AddRandomRange(const interval<double>& i) {
        random.emplace_back(i);
    }

    virtual void ManualClusterCorrection(OptionsPtr opts);
};

}} // namespace ant::expconfig
