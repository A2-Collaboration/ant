#include "Setup.h"

#include "base/Paths.h"
#include "base/Logger.h"

using namespace ant::expconfig;


std::list<std::shared_ptr<ant::Calibration::PhysicsModule> > Setup::GetCalibrations() const {
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

std::list<std::shared_ptr<ant::ReconstructHook::Base> > Setup::GetReconstructHooks() const {
    std::list< std::shared_ptr< ReconstructHook::Base > > list = reconstruct_hooks;
    for(const auto& calib : calibrations) {
        std_ext::AddToSharedPtrList<ReconstructHook::Base, Calibration::BaseModule>(
                    calib, list
                    );
    }
    for(const auto& detector : detectors) {
        std_ext::AddToSharedPtrList<ReconstructHook::Base, Detector_t>(
                    detector, list
                    );
    }
    return list;
}

std::list<std::shared_ptr<ant::Updateable_traits> > Setup::GetUpdateables() const {
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

std::string Setup::GetPIDCutsDirectory() const {
    return std::string(ANT_PATH_DATABASE)+"/"+GetName()+"/cuts";
}

std::string Setup::GetPhysicsFilesDirectory() const
{
    return std::string(ANT_PATH_DATABASE)+"/"+GetName()+"/physics_files";
}

Setup::Setup(const std::string& name, OptionsPtr opts) :
    name_(name),
    includeIgnoredElements(opts->Get<bool>("IncludeIgnoredElements", false))
{
    std::string calibrationDataFolder = std::string(ANT_PATH_DATABASE)+"/"+GetName()+"/calibration";
    calibrationDataManager = std::make_shared<calibration::DataManager>(calibrationDataFolder);
}

void Setup::BuildMappings(std::vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                          std::vector<UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
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

void Setup::IgnoreDetectorChannel(ant::Detector_t::Type_t type, unsigned channel) {
    if(includeIgnoredElements) {
        LOG_N_TIMES(1, WARNING) << "Ignored elements requested to be included nonetheless";
        return;
    }
    for(auto& detector : detectors) {
        if(detector->Type == type) {
            detector->SetIgnored(channel);
            return;
        }
    }
    LOG(WARNING) << "Setup " << GetName() << " ignored non-existing channel "
                 << channel << " of detector " << Detector_t::ToString(type);
}

void Setup::IgnoreDetectorChannels(ant::Detector_t::Type_t type, const std::vector<unsigned>& channels) {
    for(unsigned channel : channels)
        IgnoreDetectorChannel(type, channel);
}
