#pragma once

/** @page expconfig Experimental configuration
 *
 * This part represents the experimental config, including the feature to
 * automatically search and find an appropiate experimental config based on
 * ant::THeaderInfo or on the name of a setup.
 *
 */

#include "base/Detector_t.h"
#include "reconstruct/Reconstruct_traits.h"
#include "calibration/Calibration.h"
#include "base/printable.h"

#include <memory>
#include <list>
#include <map>

namespace ant {

class THeaderInfo;
class Updateable_traits;
namespace calibration {
class DataManager;
}

class ExpConfig
{
public:


    // all configs have a common base and should match via THeaderInfo
    class Base {
    public:
        virtual ~Base() = default;
        virtual bool Matches(const THeaderInfo& header) const = 0;
    };

    // the ExpConfig::Setup provides general information about the experiment
    class Setup : public virtual Base {
    public:
        static std::string ManualName;

        virtual std::string GetName() const = 0;
        virtual double GetElectronBeamEnergy() const = 0;
        virtual std::list< std::shared_ptr< Calibration::PhysicsModule> > GetCalibrations() const = 0;
        virtual std::string GetPIDCutsDirectory() const = 0;
        virtual std::shared_ptr<calibration::DataManager> GetCalibrationDataManager() const = 0;

        // you may obtain such an Expconfig::Setup via headerInfo, name,
        // get all of them, or the last found one
        static std::shared_ptr<Setup> Get(const THeaderInfo& header);
        static std::shared_ptr<Setup> Get(const std::string& name);
        static std::list<std::string> GetNames();
        static std::shared_ptr<Setup> GetLastFound();
        static std::shared_ptr<Detector_t> GetDetector(Detector_t::Type_t type);

        template<typename DetectorType>
        static std::shared_ptr<DetectorType> GetDetector();

        static void Cleanup();
    private:
        friend class ExpConfig;
        static std::shared_ptr<Setup> lastFound;
    };

    // in order to run the Reconstruction,
    // the following methods are needed
    class Reconstruct : public virtual Base {
    public:
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

        virtual candidatebuilder_config_t GetCandidateBuilderConfig() const { return candidatebuilder_config_t(); }

        // factory method to obtain such a type of config
        static std::shared_ptr<Reconstruct> Get(const THeaderInfo& header);
    };

    // each unpacker has its own config,
    // we enforce this by templates
    template<class T>
    class Unpacker : public virtual Base {
    public:
        static std::shared_ptr<T> Get(const THeaderInfo& header);
    };

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };
    class ExceptionNoConfig : public Exception {
        using Exception::Exception;
    };

    ExpConfig() = delete; // this class is more a wrapper for handling the config

private:
    template<typename T>
    static std::shared_ptr<T> Get_(const THeaderInfo& header);
};

template<typename DetectorType>
std::shared_ptr<DetectorType> ExpConfig::Setup::GetDetector()
{
    auto config = std::dynamic_pointer_cast<Reconstruct, Setup>(GetLastFound());
    if(config == nullptr)
        return nullptr;
    for(const auto& detector : config->GetDetectors()) {
        auto detector_ = std::dynamic_pointer_cast<DetectorType, Detector_t>(detector);
        if(detector_ != nullptr)
            return detector_;
    }
    return nullptr;
}

} // namespace ant
