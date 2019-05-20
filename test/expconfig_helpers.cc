#include "catch.hpp"
#include "expconfig_helpers.h"

#include "expconfig/setups/SetupRegistry.h"
#include "expconfig/ExpConfig.h"
#include "expconfig/setups/Setup.h"

#include "base/std_ext/math.h"

#include "detectors/Trigger.h"
#include "detectors/CB.h"
#include "detectors/PID.h"
#include "detectors/TAPS.h"
#include "detectors/TAPSVeto.h"
#include "detectors/EPT.h"

#include "calibration/modules/Time.h"
#include "calibration/modules/CB_Energy.h"
#include "calibration/modules/CB_TimeWalk.h"
#include "calibration/modules/PID_Energy.h"
#include "calibration/modules/PID_PhiAngle.h"
#include "calibration/modules/TAPS_Time.h"
#include "calibration/modules/TAPS_Energy.h"
#include "calibration/modules/TAPS_ShortEnergy.h"
#include "calibration/modules/TAPS_ShowerCorrection.h"
#include "calibration/modules/TAPSVeto_Energy.h"
#include "calibration/modules/TAPSVeto_Time.h"


#include "calibration/fitfunctions/FitGaus.h"
#include "calibration/fitfunctions/FitGausPol0.h"
#include "calibration/fitfunctions/FitGausPol3.h"

#include "calibration/converters/MultiHit.h"
#include "calibration/converters/MultiHitReference.h"
#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/CATCH_TDC.h"

#include "base/tmpfile_t.h"

#include <limits>

using namespace std;
using namespace ant::test;
using namespace ant::expconfig;

void ant::test::EnsureSetup(bool includeIgnored) {

    struct Setup_Test : expconfig::Setup {

        static OptionsPtr MakeOptions(bool includeIgnored) {
            auto options = make_shared<OptionsList>();
            if(includeIgnored)
                options->SetOption("IncludeIgnoredElements=1");
            return options;
        }

        tmpfolder_t calibrationDataFolder;

        Setup_Test(bool includeIgnored) : Setup("Setup_Test", MakeOptions(includeIgnored)) {

            // make sure we don't write stuff into the official database
            calibrationDataManager = std::make_shared<calibration::DataManager>(calibrationDataFolder.foldername);

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
            const bool pizzaInstalled = false;
            auto taps = make_shared<detector::TAPS_2013_11>(cherenkovInstalled, pizzaInstalled, false); // no Cherenkov, no Pizza, don't use sensitive channels
            AddDetector(taps);
            auto tapsVeto = make_shared<detector::TAPSVeto_2014>(cherenkovInstalled, pizzaInstalled); // no Cherenkov, no Pizza
            AddDetector(tapsVeto);


            // then calibrations need some rawvalues to "physical" values converters
            // they can be quite different (especially for the COMPASS TCS system), but most of them simply decode the bytes
            // to 16bit signed values
            /// \todo check if 16bit signed is correct for all those detectors
            const auto& convert_MultiHit16bit = make_shared<calibration::converter::MultiHit<std::uint16_t>>();
            const auto& convert_CATCH_Tagger = make_shared<calibration::converter::CATCH_TDC>(
                                                   trigger->Reference_CATCH_TaggerCrate
                                                   );
            const auto& convert_CATCH_CB = make_shared<calibration::converter::CATCH_TDC>(
                                               trigger->Reference_CATCH_CBCrate
                                               );
            const auto& convert_GeSiCa_SADC = make_shared<calibration::converter::GeSiCa_SADC>();
            const auto& convert_V1190_TAPSPbWO4 =  make_shared<calibration::converter::MultiHitReference<std::uint16_t>>(
                                                       trigger->Reference_V1190_TAPSPbWO4,
                                                       calibration::converter::Gains::V1190_TDC
                                                       );

            // the order of the reconstruct hooks is important
            // add both CATCH converters and the V1190 first,
            // since they need to scan the detector read for their reference hit
            AddHook(convert_CATCH_Tagger);
            AddHook(convert_CATCH_CB);
            AddHook(convert_V1190_TAPSPbWO4);

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
                                                   convert_V1190_TAPSPbWO4,
                                                   interval<double>{-500, 500}, // for BaF2
                                                   interval<double>{-500, 500}  // for PbWO4
                                                   );
            AddCalibration<calibration::TAPSVeto_Time>(tapsVeto,
                                                       calibrationDataManager,
                                                       convert_MultiHit16bit,   // for BaF2
                                                       convert_V1190_TAPSPbWO4, // for PbWO4
                                                       interval<double>{-1000, 1000},
                                                       interval<double>{-1000, 1000}
                                                       );

            AddCalibration<calibration::CB_Energy>(cb, calibrationDataManager, convert_GeSiCa_SADC,
                                                   std::vector<double>{0.0},    // default pedestal
                                                   std::vector<double>{0.07}, // default gain
                                                   std::vector<double>{2.0},    // default MeV threshold
                                                   std::vector<double>{1.0}   // default relative gain
                                                   );

            AddCalibration<calibration::PID_Energy>(pid, calibrationDataManager, convert_MultiHit16bit,
                                                    std::vector<double>{100.0},   // default pedestals
                                                    std::vector<double>{0.014},   // default gain
                                                    std::vector<double>{15.0},    // default Raw threshold
                                                    std::vector<double>{0.1},     // default MC MeV threshold
                                                    std::vector<double>{1.0}      // default relative gain
                                                    );

            AddCalibration<calibration::TAPS_Energy>(taps, calibrationDataManager, convert_MultiHit16bit,
                                                     std::vector<double>{100.0}, // default pedestal
                                                     std::vector<double>{0.3}, // default gain
                                                     5.0, 0, // default Raw thresholds BaF2/PbWO4
                                                     std::vector<double>{1.0}, // default MC MeV threshold
                                                     std::vector<double>{1.0}  // default relative gain
                                                     );

            AddCalibration<calibration::TAPS_ShortEnergy>(taps, calibrationDataManager, convert_MultiHit16bit );

            AddCalibration<calibration::TAPSVeto_Energy>(tapsVeto, calibrationDataManager, convert_MultiHit16bit);

            // enable TAPS shower correction, which is a hook running on list of clusters
            AddHook<calibration::TAPS_ShowerCorrection>();

            // the PID calibration is a physics module only
            AddCalibration<calibration::PID_PhiAngle>(pid, calibrationDataManager);

            // CB timing needs timewalk correction
            AddCalibration<calibration::CB_TimeWalk>(cb,
                                                     calibrationDataManager,
                                                     interval<double>{-std_ext::inf, std_ext::inf});
        }

        double GetElectronBeamEnergy() const override {
            return 1604.0;
        }
    };

    auto setup = make_shared<Setup_Test>(includeIgnored);
    expconfig::SetupRegistry::AddSetup(setup->GetName(), setup);
    ExpConfig::Setup::SetByName(setup->GetName());

    REQUIRE_NOTHROW(ExpConfig::Setup::Get());
}

