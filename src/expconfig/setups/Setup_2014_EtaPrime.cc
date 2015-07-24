#include "Setup.h"

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_EtaPrime :
        public Setup
{
public:
    virtual std::string GetName() const override {
        return "Setup_2014_EtaPrime";
    }

    Setup_2014_EtaPrime() {
        const auto trigger = std::make_shared<detector::Trigger>();

        const bool cherenkovInstalled = false;

        AddDetector(trigger);
        AddDetector<detector::EPT_2014>(GetBeamEnergy());
        AddDetector<detector::CB>();
        AddDetector<detector::PID_2014>();
        AddDetector<detector::TAPS_2013>(cherenkovInstalled, false); // no Cherenkov, don't use sensitive channels
        AddDetector<detector::TAPSVeto_2014>(cherenkovInstalled); // no Cherenkov

        // the order of the calibrations is important
        // since they may depend on each other
        AddCalibration<calibration::TimingCATCH>(Detector_t::Type_t::EPT,
                                                 trigger->Reference_CATCH_TaggerCrate);
        AddCalibration<calibration::TimingCATCH>(Detector_t::Type_t::CB,
                                                 trigger->Reference_CATCH_CBCrate);
        AddCalibration<calibration::IntegralSADC>(Detector_t::Type_t::CB);
        AddCalibration<calibration::TimingCATCH>(Detector_t::Type_t::PID,
                                                 trigger->Reference_CATCH_CBCrate);
        AddCalibration<calibration::IntegralCAEN>(Detector_t::Type_t::PID);
        AddCalibration<calibration::TimingTAPS>();
        AddCalibration<calibration::IntegralTAPS>();
    }

    virtual double GetBeamEnergy() const override {
        return 1604.0;
    }

    virtual cluster_thresholds_t GetClusterThresholds() const override {
        return {
            {Detector_t::Type_t::CB,   15}, // in MeV
            {Detector_t::Type_t::TAPS, 20}, // in MeV
        };
    }

    bool Matches(const THeaderInfo& header) const override {
        if(!Setup::Matches(header))
            return false;
        /// \todo Make beamtime match stricter than just detectors
        return true;

    }

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const
    {
        Setup::BuildMappings(hit_mappings, scaler_mappings);
        // you may tweak the mapping at this location here
        // for example, ignore elements
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_EtaPrime)

}}} // namespace ant::expconfig::setup
