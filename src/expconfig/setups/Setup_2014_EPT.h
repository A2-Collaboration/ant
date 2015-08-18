#include "Setup.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_EPT : public Setup
{
public:

    Setup_2014_EPT(const string& name) :
        Setup(name)
    {

        // setup the detectors of interest
        auto trigger = make_shared<detector::Trigger_2014>();
        AddDetector(trigger);

        auto EPT = make_shared<detector::EPT_2014>(GetElectronBeamEnergy());;
        AddDetector(EPT);

        auto cb = make_shared<detector::CB>();
        AddDetector(cb);

        auto pid = make_shared<detector::PID_2014>();
        AddDetector(pid);

        const bool cherenkovInstalled = false;
        auto taps = make_shared<detector::TAPS_2013>(cherenkovInstalled, false); // no Cherenkov, don't use sensitive channels
        AddDetector(taps);
        auto tapsVeto = make_shared<detector::TAPSVeto_2014>(cherenkovInstalled); // no Cherenkov
        AddDetector(tapsVeto);


        // then calibrations need some rawvalues to "physical" values converters
        // they can be quite different (especially for the COMPASS TCS system), but most of them simply decode the bytes
        // to 16bit signed values
        /// \todo check if 16bit signed is correct for all those detectors
        const auto& convert_MultiHit16bit = make_shared<calibration::converter::MultiHit16bit>();
        const auto& convert_CATCH_Tagger = make_shared<calibration::converter::CATCH_TDC>(trigger->Reference_CATCH_TaggerCrate);
        const auto& convert_CATCH_CB = make_shared<calibration::converter::CATCH_TDC>(trigger->Reference_CATCH_CBCrate);
        const auto& convert_GeSiCa_SADC = make_shared<calibration::converter::GeSiCa_SADC>();
        const auto& convert_ScalerFrequency_Beampolmon
                = make_shared<calibration::converter::ScalerFrequency>(trigger->Scaler_Beampolmon_1MHz);

        // the order of the reconstruct hooks is important
        // add both CATCH converters first,
        // since they need to scan the detector read for their reference hit
        AddHook(convert_CATCH_Tagger);
        AddHook(convert_CATCH_CB);
        // also the ScalerFrequency needs a reference
        AddHook(convert_ScalerFrequency_Beampolmon);

        AddCalibration<calibration::Scaler>(Detector_t::Type_t::EPT,
                                            convert_ScalerFrequency_Beampolmon);

        // then we add the others, and link it to the converters
        AddCalibration<calibration::Time>(EPT,
                                          calibrationDataManager,
                                          convert_CATCH_Tagger,
                                          -325 // default offset in ns
                                          );
        AddCalibration<calibration::Time>(cb,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,      // default offset in ns
                                          interval<double>{-100, 100} // default time window cut in ns
                                          );
        AddCalibration<calibration::Time>(pid,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,
                                          interval<double>{-500, 500} // default time window cut in ns
                                          );
        AddCalibration<calibration::Time>(taps,
                                          calibrationDataManager,
                                          convert_MultiHit16bit,
                                          -300, /// \todo different default for PbWO
                                          interval<double>{-500, 500},
                                          -0.100 /// \todo give measured time gains for BaF2
                                          );
        AddCalibration<calibration::Time>(tapsVeto,
                                          calibrationDataManager,
                                          convert_MultiHit16bit,
                                          160,
                                          interval<double>{-1000, 1000}, /// \todo make this window smaller...
                                          -0.05 // default gain
                                          );

        AddCalibration<calibration::CB_Energy>(cb, calibrationDataManager, convert_GeSiCa_SADC );

        AddCalibration<calibration::PID_Energy>(pid, calibrationDataManager, convert_MultiHit16bit );

        AddCalibration<calibration::TAPS_Energy>(taps, calibrationDataManager, convert_MultiHit16bit );

        AddCalibration<calibration::TAPSVeto_Energy>(calibrationDataManager, convert_MultiHit16bit);

        // enable TAPS shower correction, which is a hook running on list of clusters
        AddCalibration<calibration::TAPS_ShowerCorrection>();

        // the PID calibration is a physics module only
        AddCalibration<calibration::PID_PhiAngle>(pid, calibrationDataManager);

    }

    virtual double GetElectronBeamEnergy() const override {
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
        return true;
    }

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const
    {
        // build the mappings from the given detectors
        // that should provide sane and correct defaults
        Setup::BuildMappings(hit_mappings, scaler_mappings);

        // now you may tweak the mapping at this location here
        // for example, ignore elements
    }
};

}}} // namespace ant::expconfig::setup
