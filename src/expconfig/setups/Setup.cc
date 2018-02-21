#include "Setup.h"

#include "base/Paths.h"
#include "base/Logger.h"

#include "calibration/modules/ClusterCorrection.h"

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
    LOG_IF(includeIgnoredElements, WARNING) << "Including ignored detector elements";
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


/**
 * @brief Check options for manual modifications of cluster energies
 * searches for the string manualClusterFactor<Detector>_<Type> or manualClusterOffset<Detector>_<Type>
 * and the corresponding double value, where <Type> is the data type this value should
 * be applied to and can be MC, Data, or Both; and <Detector> is can be CB or TAPS, or can be empty, in which
 * case it will be applied to both CB and TAPS. The flags could look like
 * e.g. manualClusterFactorCB_MC=1.1 or manualClusterOffsetTAPS_Both=5 or manualClusterOffset_Data=4.5
 *
 * @param OptionsPtr
 * @see ClusterCorrFactor and ClusterCorrOffset
 */
void Setup::ManualClusterCorrection(OptionsPtr opts)
{
    if (!opts->HasOptionStartsWith("manualCluster"))
        return;

    using Filter_t = calibration::ClusterCorrection::Filter_t;

    std::string apply_str;

    auto get_filter = [&apply_str] (const std::string opt) -> Filter_t {
        if (std_ext::contains(opt, "_MC")) {
            apply_str = "MC";
            return Filter_t::MC;
        } else if (std_ext::contains(opt, "_Data")) {
            apply_str = "Data";
            return Filter_t::Data;
        } else if (std_ext::contains(opt, "_Both")) {
            apply_str = "both Data and MC";
            return Filter_t::Both;
        } else {
            LOG(WARNING) << "No filter given for option '" << opt << "'; assume both Data and MC";
            apply_str = "both Data and MC";
            return Filter_t::Both;
        }
    };

    auto get_detector = [this] (const Detector_t::Type_t det_type) -> std::shared_ptr<Detector_t> {
        for (const auto& detector : GetDetectors()) {
            if (detector->Type == det_type)
                return detector;
        }
        throw ExpConfig::ExceptionNoDetector("Could not find detector");
    };

    std::shared_ptr<ClusterDetector_t> det;

    auto determine_detectors = [&det, &get_detector] (const std::string opt) {
        if (std_ext::contains(opt, "CB"))
            det = std::dynamic_pointer_cast<ClusterDetector_t>(get_detector(Detector_t::Type_t::CB));
        else if (std_ext::contains(opt, "TAPS"))
            det = std::dynamic_pointer_cast<ClusterDetector_t>(get_detector(Detector_t::Type_t::TAPS));
        else {
            LOG(WARNING) << "No detector given in option '" << opt << "'; assume both CB and TAPS";
            return true;
        }

        return false;
    };

    std::string option;
    bool CBandTAPS;

    while (opts->HasUnusedOptionStartsWith("manualCluster")) {
        CBandTAPS = false;
        // check for cluster energy correction factor
        option = opts->UnusedOptionStartsWith("manualClusterFactor");
        if (!option.empty()) {
            const double value = opts->Get<double>(option, 1.);
            auto filter = get_filter(option);

            CBandTAPS = determine_detectors(option);

            if (!CBandTAPS) {
                LOG(INFO) << "Apply cluster correction factor of " << value << " to "
                          << Detector_t::ToString(det->Type) << " on " << apply_str;
                AddCalibration<calibration::ClusterCorrFactor>(det, "ManualClusterEFactor",
                                                               filter, nullptr, value);
            } else {
                LOG(INFO) << "Apply cluster correction factor of " << value
                          << " to both CB and TAPS on " << apply_str;
                for (const auto d : {"CB", "TAPS"}) {
                    determine_detectors(d);
                    AddCalibration<calibration::ClusterCorrFactor>(det, "ManualClusterEFactor",
                                                                   filter, nullptr, value);
                }
            }
        }

        // check for cluster energy offset correction
        option = opts->UnusedOptionStartsWith("manualClusterOffset");
        if (!option.empty()) {
            const double value = opts->Get<double>(option, 0.);
            auto filter = get_filter(option);

            CBandTAPS = determine_detectors(option);

            if (!CBandTAPS) {
                LOG(INFO) << "Apply cluster offset correction of " << value << " MeV to "
                          << Detector_t::ToString(det->Type) << " on " << apply_str;
                AddCalibration<calibration::ClusterCorrOffset>(det, "ManualClusterEOffset",
                                                               filter, nullptr, value);
            } else {
                LOG(INFO) << "Apply cluster offset correction of " << value
                          << " MeV to both CB and TAPS on " << apply_str;
                for (const auto d : {"CB", "TAPS"}) {
                    determine_detectors(d);
                    AddCalibration<calibration::ClusterCorrOffset>(det, "ManualClusterEOffset",
                                                                   filter, nullptr, value);
                }
            }
        }
    }

}
