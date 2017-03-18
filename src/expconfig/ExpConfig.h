#pragma once

/** @page expconfig Experimental configuration
 *
 * This part represents the experimental config, including the feature to
 * automatically search and find an appropiate experimental config based on
 * ant::TID or on the name of a setup.
 *
 */

#include "base/Detector_t.h"
#include "reconstruct/Reconstruct_traits.h"
#include "calibration/Calibration.h"
#include "base/printable.h"
#include "base/piecewise_interval.h"

#include <memory>
#include <list>
#include <map>

namespace ant {

struct TID;
struct Updateable_traits;
namespace calibration {
class DataManager;
}

namespace expconfig {
class Setup_traits;
} // keep ExpConfig itself in ant namespace

class ExpConfig
{
public:


    using SetupPtr = std::shared_ptr<expconfig::Setup_traits>;

    class Setup {
    public:

        /**
         * @brief Get
         * @return
         */
        static SetupPtr Get();

        template<typename SetupType = expconfig::Setup_traits>
        static std::shared_ptr<SetupType> GetByType();

        /**
         * @brief GetDetector ask for detector by type
         * @param type the requested type
         * @return detector pointer, never nullptr
         * @throws ExceptionNoDetector if no detector was found
         * @see the templated version to get type-safe detector
         */
        static std::shared_ptr<Detector_t> GetDetector(Detector_t::Type_t type);

        template<typename DetectorType>
        static std::shared_ptr<DetectorType> GetDetector();

        // give this registry an idea how to search for a setup
        static void SetByName(const std::string& setupname);
        static void SetByTID(const TID& tid);

        static std::list<std::string> GetNames();
        static void Cleanup();

    private:
        friend class ExpConfig;
        static std::string manualName;
        static ExpConfig::SetupPtr currentSetup;
    };

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };
    class ExceptionNoSetup : public Exception {
        using Exception::Exception;
    };
    class ExceptionNoDetector : public Exception {
        using Exception::Exception;
    };
    ExpConfig() = delete; // this class is more a wrapper for handling the config

};

namespace expconfig {

/**
 * @brief The Setup_traits class is the interface to the "static" experimental information
 */
class Setup_traits {
public:
    virtual bool Matches(const TID& header) const = 0;

    virtual std::string GetName() const = 0;
    virtual double GetElectronBeamEnergy() const = 0;
    virtual std::list< std::shared_ptr< Calibration::PhysicsModule> > GetCalibrations() const = 0;
    virtual std::string GetPIDCutsDirectory() const = 0;
    virtual std::string GetPhysicsFilesDirectory() const = 0;
    virtual std::shared_ptr<calibration::DataManager> GetCalibrationDataManager() const = 0;

    virtual std::list< std::shared_ptr< ReconstructHook::Base > > GetReconstructHooks() const = 0;
    virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const = 0;
    virtual std::list< std::shared_ptr< Updateable_traits> > GetUpdateables() const = 0;

    struct candidatebuilder_config_t {
        /// @see ant::reconstruct::CandidateBuilder::PID_phi_epsilon
        double PID_Phi_Epsilon = 0.0;      // in Rad
        double CB_ClusterThreshold = 15;   // in MeV
        double TAPS_ClusterThreshold = 20; // in MeV
        candidatebuilder_config_t() = default;
    };

    virtual candidatebuilder_config_t GetCandidateBuilderConfig() const {
        // return defaults by default :)
        return candidatebuilder_config_t();
    }

    virtual bool GetIncludeIgnoredElements() const = 0;
    virtual ant::PiecewiseInterval<double> GetPromptWindows() const = 0;
    virtual ant::PiecewiseInterval<double> GetRandomWindows() const = 0;

    virtual ~Setup_traits() = default;

}; // Setup_traits

} // namespace ant::expconfig


// templated static getters can now be implemented,
// as Setup_traits is now defined

template<typename DetectorType>
std::shared_ptr<DetectorType> ExpConfig::Setup::GetDetector()
{
    auto setup = Get();
    if(setup == nullptr)
        throw ExceptionNoSetup("Could not find setup to search for required detector");
    for(const auto& detector : setup->GetDetectors()) {
        auto detector_ = std::dynamic_pointer_cast<DetectorType, Detector_t>(detector);
        if(detector_ != nullptr)
            return detector_;
    }
    throw ExceptionNoDetector("Could not find detector in given setup");
}

template<typename SetupType>
std::shared_ptr<SetupType> ExpConfig::Setup::GetByType()
{
    auto setup = std::dynamic_pointer_cast<SetupType, expconfig::Setup_traits>(currentSetup);
    if(!setup)
        throw ExceptionNoSetup("No setup found with requested type");
    return setup;
}



} // namespace ant
