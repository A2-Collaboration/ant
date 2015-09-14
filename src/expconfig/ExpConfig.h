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

class ExpConfig
{
public:

    static std::string ManualSetupName;

    // all configs have a common base and should match via THeaderInfo
    class Base {
    public:
        virtual ~Base() = default;
        virtual bool Matches(const THeaderInfo& header) const = 0;
    };

    // the ExpConfig::Module provides general information about the experiment
    class Setup : public virtual Base {
    public:
        virtual std::string GetName() const = 0;
        virtual double GetElectronBeamEnergy() const = 0;
        virtual std::list< std::shared_ptr< Calibration::PhysicsModule> > GetCalibrations() const = 0;
        virtual std::string GetPIDCutsDirectory() const = 0;

        // you may obtain such an Expconfig::Module via headerInfo, name,
        // get all of them, or the last found one
        static std::shared_ptr<Setup> Get(const THeaderInfo& header);
        static std::shared_ptr<Setup> Get(const std::string& name);
        static std::list<std::string> GetNames();
        static std::shared_ptr<Setup> GetLastFound();

        static void Cleanup();
    private:
        static std::list< std::shared_ptr<Setup> > getAll();
    };

    // in order to run the Reconstruction,
    // the following methods are needed
    class Reconstruct : public virtual Base {
    public:
        virtual std::list< std::shared_ptr< ReconstructHook::Base > > GetReconstructHooks() const = 0;
        virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const = 0;
        virtual std::list< std::shared_ptr< Updateable_traits> > GetUpdateables() const = 0;

        // for clustering, may be extended to config struct
        // (e.g. for position weighting options)
        using cluster_thresholds_t = std::map<Detector_t::Type_t, double>;
        virtual cluster_thresholds_t  GetClusterThresholds() const = 0;

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
    static std::shared_ptr<Setup> lastSetupFound;

};

} // namespace ant
