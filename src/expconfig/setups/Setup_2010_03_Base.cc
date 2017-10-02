#include "Setup.h"

#include "Setup_2010_03_Base.h"

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

#include "TH1.h"
#include "TF1.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {


Setup_2010_03_Base::Setup_2010_03_Base(const std::string& name, OptionsPtr opt) : Setup(name, opt),
    MCTaggerHits(opt->Get<bool>("MCTaggerHits",false)),
    cherenkovInstalled(false),
    Trigger(make_shared<detector::Trigger_2007>()),
    Tagger(make_shared<detector::Tagger_2010_03>()),
    CB(make_shared<detector::CB>()),
    PID(make_shared<detector::PID_2009_07>()),
    TAPS(make_shared<detector::TAPS_2009_03>(cherenkovInstalled, false)), // no cherenkov, don't use sensitive channels
    TAPSVeto(make_shared<detector::TAPSVeto_2009_03>(cherenkovInstalled)) // no cherenkov
{
        AddDetector(CB);
        AddDetector(PID);
        AddDetector(Trigger);
        AddDetector(Tagger);
        AddDetector(TAPS);
        AddDetector(TAPSVeto);

        // then calibrations need some rawvalues to "physical" values converters
        // they can be quite different (especially for the COMPASS TCS system), but most of them simply decode the bytes
        // to 16bit signed values
        /// \todo check if 16bit unsigned is correct for all those detectors
        const auto& convert_MultiHit16bit = make_shared<calibration::converter::MultiHit<std::uint16_t>>();
        const auto& convert_CATCH_Tagger = make_shared<calibration::converter::CATCH_TDC>(
                                               Trigger->Reference_CATCH_TaggerCrate
                                               );
        const auto& convert_CATCH_CB = make_shared<calibration::converter::CATCH_TDC>(
                                           Trigger->Reference_CATCH_CBCrate
                                           );
        const auto& convert_GeSiCa_SADC = make_shared<calibration::converter::GeSiCa_SADC>();
        const auto& convert_V1190_TAPSPbWO4 =  make_shared<calibration::converter::MultiHitReference<std::uint16_t>>(
//                                                   Trigger->Reference_V1190_TAPSPbWO4,
                                                   Trigger->Reference_CATCH_CBCrate,
                                                   calibration::converter::Gains::V1190_TDC
                                                   );

        // the order of the reconstruct hooks is important
        // add both CATCH converters and the V1190 first,
        // since they need to scan the detector read for their reference hit
        AddHook(convert_CATCH_Tagger);
        AddHook(convert_CATCH_CB);
        AddHook(convert_V1190_TAPSPbWO4);

        // Tagging efficiencies are loaded via a calibration module
//        AddCalibration<calibration::TaggEff>(Tagger, calibrationDataManager);

        const bool timecuts = !opt->Get<bool>("DisableTimecuts");
        interval<double> no_timecut(-std_ext::inf, std_ext::inf);
        if(!timecuts)
            LOG(INFO) << "Disabling timecuts";

        const bool thresholds = !opt->Get<bool>("DisableThresholds");
        if(!thresholds)
            LOG(INFO) << "Disabling thresholds";

        // then we add the others, and link it to the converters
        AddCalibration<calibration::Time>(Tagger,
                                          calibrationDataManager,
                                          convert_CATCH_Tagger,
                                          -325, // default offset in ns
                                          std::make_shared<calibration::gui::FitGausPol0>(),
                                          timecuts ? interval<double>{-120, 120} : no_timecut
                                          );
        AddCalibration<calibration::Time>(CB,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,      // default offset in ns
                                          std::make_shared<calibration::gui::CBPeakFunction>(),
                                          // Let CB_TimeWalk decide on good timing hits
                                          // as there are some broken TDCs which may be recovered
                                          // from energy information
                                          no_timecut
                                          );
        AddCalibration<calibration::Time>(PID,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,
                                          std::make_shared<calibration::gui::FitGaus>(),
                                          // The PID timing must be plotted on a "clean" sample versus
                                          // energy, for example identify good pi0 events with protons in CB
                                          // with kinematic fitter. See ProtonPi0 physics class.
                                          timecuts ? interval<double>{-25, 40} : no_timecut
                                          );
        AddCalibration<calibration::TAPS_Time>(TAPS,
                                               calibrationDataManager,
                                               convert_MultiHit16bit,   // for BaF2
                                               convert_V1190_TAPSPbWO4, // for PbWO4
                                               timecuts ? interval<double>{-15, 15} : no_timecut, // for BaF2
                                               timecuts ? interval<double>{-25, 25} : no_timecut  // for PbWO4
                                               );
        AddCalibration<calibration::TAPSVeto_Time>(TAPSVeto,
                                                   calibrationDataManager,
                                                   convert_MultiHit16bit,   // for BaF2
                                                   convert_V1190_TAPSPbWO4, // for PbWO4
                                                   timecuts ? interval<double>{-12, 12} : no_timecut,
                                                   timecuts ? interval<double>{-12, 12} : no_timecut
                                                   );

        AddCalibration<calibration::CB_Energy>(CB, calibrationDataManager, convert_GeSiCa_SADC,
                                               std::vector<double>{0},    // default pedestal
                                               std::vector<double>{0.07}, // default gain
                                               // default threshold, only used on MC
                                               std::vector<double>{thresholds ? 1.2 : 0.0},
                                               std::vector<double>{1.0}   // default relative gain
                                               );

        AddCalibration<calibration::PID_Energy>(PID, calibrationDataManager, convert_MultiHit16bit,
                                                std::vector<double>{100.0},   // default pedestals
                                                std::vector<double>{0.014},   // default gain
                                                std::vector<double>{thresholds ? 15.0 : -std_ext::inf}, // default Raw threshold
                                                std::vector<double>{thresholds ? 0.1  : 0.0},     // default MC MeV threshold
                                                std::vector<double>{1.0}      // default relative gain
                                                );

        AddCalibration<calibration::TAPS_Energy>(TAPS, calibrationDataManager, convert_MultiHit16bit,
                                                 std::vector<double>{100}, // default pedestal
                                                 std::vector<double>{0.3}, // default gain
                                                 (thresholds ? 5.0 : -std_ext::inf), 0, // default Raw thresholds BaF2/PbWO4
                                                 std::vector<double>{thresholds ? 3.4 : 0.0}, // default MC MeV threshold
                                                 std::vector<double>{1.0}  // default relative gain
                                                 );

        AddCalibration<calibration::TAPS_ShortEnergy>(TAPS, calibrationDataManager, convert_MultiHit16bit );

        AddCalibration<calibration::TAPSVeto_Energy>(TAPSVeto, calibrationDataManager, convert_MultiHit16bit);

        // enable TAPS shower correction, which is a hook running on list of clusters
        AddHook<calibration::TAPS_ShowerCorrection>();

        // add ToF timing to TAPS clusters
        AddCalibration<calibration::TAPS_ToF>(TAPS, calibrationDataManager);

        // the PID calibration is a physics module only
        AddCalibration<calibration::PID_PhiAngle>(PID, calibrationDataManager);

        // CB timing needs timewalk correction
        AddCalibration<calibration::CB_TimeWalk>(CB, calibrationDataManager,
                                                 timecuts ? interval<double>{-25, 25} : no_timecut,
                                                 7 // energy threshold for BadTDCs
                                                 );



        //Cluster Smearing, Energy. Only activates if root file with histogram present in calibration data folder.
        //Place a file in the MC folder to use MC smearing. Do not put one in the "Data" calibration folder unless
        //you want to smear data as well (probably not...)

        // MC scaling was found to be superfluous, after using "clean" clusters not touching any hole
//        AddCalibration<calibration::ClusterSmearing>(CB,   "ClusterSmearing",  calibration::ClusterCorrection::Filter_t::MC, calibrationDataManager);
//        AddCalibration<calibration::ClusterSmearing>(TAPS, "ClusterSmearing",  calibration::ClusterCorrection::Filter_t::MC, calibrationDataManager);

        // ECorr, should be applied after MC smearing
//        AddCalibration<calibration::ClusterECorr>(CB,   "ClusterECorr",  calibration::ClusterCorrection::Filter_t::Both, calibrationDataManager);
//        AddCalibration<calibration::ClusterECorr>(TAPS, "ClusterECorr",  calibration::ClusterCorrection::Filter_t::Both, calibrationDataManager);

        // prompt is chosen with TriggerSimulation::GetCorrectedTaggerTime
        AddPromptRange({  -3,   3});
        AddRandomRange({ -50,  -5});
        AddRandomRange({  5,   50});
    }

double Setup_2010_03_Base::GetElectronBeamEnergy() const {
        return 450.0;
}

Setup_traits::candidatebuilder_config_t Setup_2010_03_Base::GetCandidateBuilderConfig() const {
    candidatebuilder_config_t conf;
    conf.PID_Phi_Epsilon = std_ext::degree_to_radian(2.0);
    conf.CB_ClusterThreshold = 15;
    conf.TAPS_ClusterThreshold = 20;
    return conf;
}


//AUTO_REGISTER_SETUP(Setup_2010_03_Base)


}}} // namespace ant::expconfig::setup
