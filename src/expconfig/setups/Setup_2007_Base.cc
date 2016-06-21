#include "Setup.h"

#include "detectors/CB.h"
#include "detectors/PID.h"
#include "detectors/Tagger.h"
#include "detectors/Trigger.h"
#include "detectors/TAPS.h"
#include "detectors/TAPSVeto.h"

#include "calibration/modules/Time.h"
#include "calibration/modules/CB_Energy.h"
#include "calibration/modules/CB_TimeWalk.h"
#include "calibration/modules/PID_Energy.h"
#include "calibration/modules/PID_PhiAngle.h"
#include "calibration/modules/TAPS_Time.h"
#include "calibration/modules/TAPS_Energy.h"
#include "calibration/modules/TAPS_ShortEnergy.h"
#include "calibration/modules/TAPS_ShowerCorrection.h"
#include "calibration/modules/TAPS_ToF.h"
#include "calibration/modules/TAPSVeto_Energy.h"
#include "calibration/modules/TAPSVeto_Time.h"


#include "calibration/fitfunctions/FitGaus.h"
#include "calibration/fitfunctions/FitGausPol0.h"
#include "calibration/fitfunctions/FitGausPol3.h"

#include "calibration/converters/MultiHit.h"
#include "calibration/converters/MultiHitReference.h"
#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/CATCH_TDC.h"

#include "base/std_ext/math.h"
#include "base/Logger.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2007_Base : public Setup
{
public:

    Setup_2007_Base(const std::string& name, OptionsPtr opt) : Setup(name, opt)
    {
        auto cb = make_shared<detector::CB>();
        AddDetector(cb);

        auto pid = make_shared<detector::PID_2007>();
        AddDetector(pid);

        auto trigger = make_shared<detector::Trigger_2007>();
        AddDetector(trigger);

        auto tagger = make_shared<detector::Tagger_2007>();
        AddDetector(tagger);

//        auto taps = make_shared<detector::TAPS_2007>(false, false); // no Cherenkov, don't use sensitive channels
//        AddDetector(taps);

        auto tapsVeto = make_shared<detector::TAPSVeto_2007>(false); // no Cherenkov
        AddDetector(tapsVeto);

        // then calibrations need some rawvalues to "physical" values converters
        // they can be quite different (especially for the COMPASS TCS system), but most of them simply decode the bytes
        // to 16bit signed values
        /// \todo check if 16bit unsigned is correct for all those detectors
        const auto& convert_MultiHit16bit = make_shared<calibration::converter::MultiHit<std::uint16_t>>();
        const auto& convert_CATCH_Tagger = make_shared<calibration::converter::CATCH_TDC>(
                                               trigger->Reference_CATCH_TaggerCrate
                                               );
        const auto& convert_CATCH_CB = make_shared<calibration::converter::CATCH_TDC>(
                                           trigger->Reference_CATCH_CBCrate
                                           );
        const auto& convert_GeSiCa_SADC = make_shared<calibration::converter::GeSiCa_SADC>();

        // the order of the reconstruct hooks is important
        // add both CATCH converters and the V1190 first,
        // since they need to scan the detector read for their reference hit
        AddHook(convert_CATCH_Tagger);
        AddHook(convert_CATCH_CB);


        const bool timecuts = !opt->Get<bool>("DisableTimecuts");
        interval<double> no_timecut(-std_ext::inf, std_ext::inf);
        if(!timecuts)
            LOG(INFO) << "Disabling timecuts";

        const bool thresholds = !opt->Get<bool>("DisableThresholds");
        if(!thresholds)
            LOG(INFO) << "Disabling thresholds";

        // then we add the others, and link it to the converters
        AddCalibration<calibration::Time>(tagger,
                                          calibrationDataManager,
                                          convert_CATCH_Tagger,
                                          -325, // default offset in ns
                                          std::make_shared<calibration::gui::FitGausPol0>(),
                                          timecuts ? interval<double>{-120, 120} : no_timecut
                                          );
        AddCalibration<calibration::Time>(cb,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,      // default offset in ns
                                          std::make_shared<calibration::gui::CBPeakFunction>(),
                                          // before timewalk correction
                                          timecuts ? interval<double>{-20, 200} : no_timecut
                                          );
        AddCalibration<calibration::Time>(pid,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,
                                          std::make_shared<calibration::gui::FitGaus>(),
                                          timecuts ? interval<double>{-20, 20} : no_timecut
                                          );
        AddCalibration<calibration::TAPSVeto_Time>(tapsVeto,
                                                   calibrationDataManager,
                                                   convert_MultiHit16bit,   // for BaF2
                                                   nullptr,    // for PbWO4
                                                   timecuts ? interval<double>{-12, 12} : no_timecut
                                                   );

        AddCalibration<calibration::CB_Energy>(cb, calibrationDataManager, convert_GeSiCa_SADC,
                                               0,    // default pedestal
                                               0.07, // default gain
                                               thresholds ? 2 : 0,    // default threshold
                                               1.0   // default relative gain
                                               );

        // CB timing needs timewalk correction
        AddCalibration<calibration::CB_TimeWalk>(cb, calibrationDataManager,
                                                 timecuts ? interval<double>{-10, 20} : no_timecut);

        AddCalibration<calibration::PID_Energy>(pid, calibrationDataManager, convert_MultiHit16bit );

        // the PID calibration is a physics module only
        AddCalibration<calibration::PID_PhiAngle>(pid, calibrationDataManager);

        AddCalibration<calibration::TAPSVeto_Energy>(tapsVeto, calibrationDataManager, convert_MultiHit16bit);

    }

    virtual double GetElectronBeamEnergy() const override {
        return 1508.0;
    }

    bool Matches(const TID& tid) const override
    {
        /// \todo This is for all of 2007 for now. for testing.
        if(!std_ext::time_between(tid.Timestamp, "2007-01-016", "2007-12-31"))
            return false;
        return true;
    }


    virtual ExpConfig::Setup::candidatebuilder_config_t GetCandidateBuilderConfig() const override {
        candidatebuilder_config_t conf;
        conf.PID_Phi_Epsilon = std_ext::degree_to_radian(2.0);
        conf.CB_ClusterThreshold = 15;
        conf.TAPS_ClusterThreshold = 20;
        return conf;
    }
};

AUTO_REGISTER_SETUP(Setup_2007_Base)


}}} // namespace ant::expconfig::setup
