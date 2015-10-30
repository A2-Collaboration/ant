#include "Setup.h"

#include "base/std_ext/math.h"

#include "detectors/Trigger.h"
#include "detectors/CB.h"
#include "detectors/PID.h"
#include "detectors/TAPS.h"
#include "detectors/TAPSVeto.h"
#include "detectors/EPT.h"

#include "calibration/modules/Time.h"
#include "calibration/modules/Scaler.h"
#include "calibration/modules/CB_Energy.h"
#include "calibration/modules/CB_TimeWalk.h"
#include "calibration/modules/PID_Energy.h"
#include "calibration/modules/PID_PhiAngle.h"
#include "calibration/modules/TAPS_Time.h"
#include "calibration/modules/TAPS_Energy.h"
#include "calibration/modules/TAPS_ShortEnergy.h"
#include "calibration/modules/TAPS_ShowerCorrection.h"
#include "calibration/modules/TAPSVeto_Energy.h"

#include "calibration/fitfunctions/FitGaus.h"
#include "calibration/fitfunctions/FitGausPol0.h"
#include "calibration/fitfunctions/FitGausPol3.h"

#include "calibration/converters/CATCH_TDC.h"
#include "calibration/converters/MultiHit16bit.h"
#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/ScalerFrequency.h"

#include <limits>

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief The Setup_Test class exists to run the test cases below test/
 */
class Setup_Test :
        public Setup
{
public:

    Setup_Test(const std::string& name, SetupOptPtr opt) : Setup(name, opt) {
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
                                          -325, // default offset in ns
                                          std::make_shared<calibration::gui::FitGausPol0>()
                                          );
        AddCalibration<calibration::Time>(cb,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,      // default offset in ns
                                          std::make_shared<calibration::gui::FitGaus>(),
                                          interval<double>{-60, 60} // default time window cut in ns
                                          );
        AddCalibration<calibration::Time>(pid,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,
                                          std::make_shared<calibration::gui::FitGaus>(),
                                          interval<double>{-500, 500} // default time window cut in ns
                                          );
        AddCalibration<calibration::TAPS_Time>(taps,
                                          calibrationDataManager,
                                          convert_MultiHit16bit,
                                          interval<double>{-500, 500}
                                          );
        AddCalibration<calibration::Time>(tapsVeto,
                                          calibrationDataManager,
                                          convert_MultiHit16bit,
                                          -100,
                                          std::make_shared<calibration::gui::FitGausPol0>(),
                                          interval<double>{-1000, 1000}, /// \todo make this window smaller...
                                          -0.05 // default gain
                                          );

        AddCalibration<calibration::CB_Energy>(cb, calibrationDataManager, convert_GeSiCa_SADC );

        AddCalibration<calibration::PID_Energy>(pid, calibrationDataManager, convert_MultiHit16bit );

        AddCalibration<calibration::TAPS_Energy>(taps, calibrationDataManager, convert_MultiHit16bit );

        AddCalibration<calibration::TAPS_ShortEnergy>(taps, calibrationDataManager, convert_MultiHit16bit );

        AddCalibration<calibration::TAPSVeto_Energy>(tapsVeto, calibrationDataManager, convert_MultiHit16bit);

        // enable TAPS shower correction, which is a hook running on list of clusters
        AddCalibration<calibration::TAPS_ShowerCorrection>();

        // the PID calibration is a physics module only
        AddCalibration<calibration::PID_PhiAngle>(pid, calibrationDataManager);

        // CB timing needs timewalk correction
        AddCalibration<calibration::CB_TimeWalk>(cb,
                                                 calibrationDataManager,
                                                 interval<double>{-std_ext::inf, std_ext::inf});
    }

    virtual double GetElectronBeamEnergy() const override {
        // don't be afraid to use NaN if you don't know a value
        return numeric_limits<double>::quiet_NaN();
    }

    bool Matches(const THeaderInfo&) const override {
        // Setup must be manually selected (for test cases)
        // via ExpConfig::Setup::ManualName
        return false;
    }

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const override
    {
        // build mappings according to detectors
        Setup::BuildMappings(hit_mappings, scaler_mappings);
    }
};

// if you forget this, your setup is never found....
AUTO_REGISTER_SETUP(Setup_Test)

}}} // namespace ant::expconfig::setup
