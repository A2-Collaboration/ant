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
        auto convert_CATCH_Tagger = make_shared<calibration::converter::CATCH_TDC>(trigger->Reference_CATCH_TaggerCrate);
        auto convert_CATCH_CB = make_shared<calibration::converter::CATCH_TDC>(trigger->Reference_CATCH_CBCrate);
        auto convert_MultiHit16bit = make_shared<calibration::converter::MultiHit16bit>();


        // add both CATCH converters first,
        // since they need to scan the detector read for their reference hit
        AddCalibration(convert_CATCH_Tagger);
        AddCalibration(convert_CATCH_CB);

        // then we add the others, and link it to the CATCH converters
        AddCalibration<calibration::Timing>(Detector_t::Type_t::EPT, convert_CATCH_Tagger);
        AddCalibration<calibration::Timing>(Detector_t::Type_t::CB,  convert_CATCH_CB);
        AddCalibration<calibration::Timing>(Detector_t::Type_t::PID, convert_CATCH_CB);
        AddCalibration<calibration::Timing>(Detector_t::Type_t::TAPS, convert_MultiHit16bit);
//        AddCalibration<calibration::Timing>(Detector_t::Type_t::TAPSVeto, convert_TAPS_Timing);


        AddCalibration<calibration::Integral>(Detector_t::Type_t::CB,
                                              make_shared<calibration::converter::GeSiCa_SADC>(),
                                              0,    // default pedestal in raw
                                              0.07, // default gain
                                              2     // default threshold in MeV
                                              );

        AddCalibration<calibration::Integral>(Detector_t::Type_t::PID,
                                              convert_MultiHit16bit,
                                              100,    // default pedestal in raw
                                              0.014,  // default gain
                                              0.001   // default threshold in MeV
                                              );

        AddCalibration<calibration::Integral>(Detector_t::Type_t::TAPS,
                                              convert_MultiHit16bit,
                                              100,   // default pedestal in raw
                                              0.30,  // default gain
                                              1      // default threshold in MeV
                                              );

//        AddCalibration<calibration::Integral>(Detector_t::Type_t::TAPSVeto,
//                                              convert_MultiHit16bit,
//                                              100,     // default pedestal in raw
//                                              0.010, // default gain
//                                              0.1    // default threshold in MeV
//                                              );
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
