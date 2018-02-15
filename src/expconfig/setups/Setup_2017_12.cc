#include "Setup_2017_12.h"

#include "base/std_ext/math.h"
#include "base/Logger.h"

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
#include "calibration/modules/Tagger_QDC.h"
#include "calibration/modules/TaggEff.h"
#include "calibration/modules/NewTagger_Time.h"
#include "calibration/modules/ClusterCorrection.h"

#include "calibration/fitfunctions/FitGaus.h"
#include "calibration/fitfunctions/FitGausPol0.h"
#include "calibration/fitfunctions/FitGausPol3.h"

#include "calibration/converters/MultiHit.h"
#include "calibration/converters/MultiHitReference.h"
#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/CATCH_TDC.h"


using namespace std;
using namespace ant::expconfig;
using namespace ant::expconfig::setup;

// OBS!!! Preliminary!! E.g. the tagger doesn't work yet

Setup_2017_12::Setup_2017_12(const string& name, OptionsPtr opt) :
    Setup(name, opt),
    MCTaggerHits(opt->Get<bool>("MCTaggerHits",false)),
    cherenkovInstalled(false),
    pizzaInstalled(false),
    Trigger(make_shared<detector::Trigger_2014>()),
    Tagger(make_shared<detector::Tagger_2017_12>()),
    CB(make_shared<detector::CB>()),
    PID(make_shared<detector::PID_2014>()),
    TAPS(make_shared<detector::TAPS_2013_11>(cherenkovInstalled, pizzaInstalled, false)), // false = don't use sensitive channels
    TAPSVeto(make_shared<detector::TAPSVeto_2014>(cherenkovInstalled, pizzaInstalled))
{
    // add the detectors of interest
    AddDetector(Trigger);
    AddDetector(Tagger);
    AddDetector(CB);
    AddDetector(PID);
    AddDetector(TAPS);
    AddDetector(TAPSVeto);

    // Possible: can set set inner ring and outer ring to NoCalib for TAPS (see 2014) "touches hole"
    // removed

    // Broken, BadTDC or NoCalib elements
    //CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {});
    //CB->SetElementFlag(Detector_t::ElementFlag_t::BadTDC, {});
    //CB->SetElementFlag(Detector_t::ElementFlag_t::NoCalibFill, {});

    // then calibrations need some rawvalues to "physical" values converters
    // they can be quite different (especially for the COMPASS TCS system), but most of them simply decode the bytes
    // to 16bit signed values
    const auto& convert_MultiHit16bit = make_shared<calibration::converter::MultiHit<std::uint16_t>>();
    // converters for the new tagger, one for each reference channel
    const auto& convert_V1190_Tagger1 = make_shared<calibration::converter::MultiHitReference<std::uint16_t>>(
                                           Trigger->Reference_V1190_TaggerTDC1,
                                           calibration::converter::Gains::V1190_TDC
                                           );
    const auto& convert_V1190_Tagger2 = make_shared<calibration::converter::MultiHitReference<std::uint16_t>>(
                                           Trigger->Reference_V1190_TaggerTDC2,
                                           calibration::converter::Gains::V1190_TDC
                                           );
    const auto& convert_V1190_Tagger3 = make_shared<calibration::converter::MultiHitReference<std::uint16_t>>(
                                           Trigger->Reference_V1190_TaggerTDC3,
                                           calibration::converter::Gains::V1190_TDC
                                           );

    const auto& convert_CATCH_CB = make_shared<calibration::converter::CATCH_TDC>(
                                       Trigger->Reference_CATCH_CBCrate
                                       );
    const auto& convert_GeSiCa_SADC = make_shared<calibration::converter::GeSiCa_SADC>();
    const auto& convert_V1190_TAPSPbWO4 =  make_shared<calibration::converter::MultiHitReference<std::uint16_t>>(
                                               Trigger->Reference_V1190_TAPSPbWO4,
                                               calibration::converter::Gains::V1190_TDC
                                               );

    // the order of the reconstruct hooks is important
    // add both CATCH converters and the V1190 first,
    // since they need to scan the detector read for their reference hit
    AddHook(convert_V1190_Tagger1);
    AddHook(convert_V1190_Tagger2);
    AddHook(convert_V1190_Tagger3);
    AddHook(convert_CATCH_CB);
    AddHook(convert_V1190_TAPSPbWO4);

    // Tagger/EPT QDC measurements need some simple hook
    AddHook<calibration::Tagger_QDC>(Tagger->Type, convert_MultiHit16bit);

    // Tagging efficiencies are loaded via a calibration module
    AddCalibration<calibration::TaggEff>(Tagger, calibrationDataManager);

    const bool timecuts = !opt->Get<bool>("DisableTimecuts");
    interval<double> no_timecut(-std_ext::inf, std_ext::inf);
    if(!timecuts)
        LOG(INFO) << "Disabling timecuts";

    const bool thresholds = !opt->Get<bool>("DisableThresholds");
    if(!thresholds)
        LOG(INFO) << "Disabling thresholds";

    const bool pedestals = !opt->Get<bool>("SetPedestalsToZero");
    if(!pedestals)
        LOG(INFO) << "Setting Pedestals for PID/TAPS/TAPSShort/TAPSVeto to 0";

    // then we add the others, and link it to the converters
    AddCalibration<calibration::NewTagger_Time>(Tagger,
                                      calibrationDataManager,
                                      std::map<detector::Tagger::TDCSector_t, Calibration::Converter::ptr_t>{
                                        {detector::Tagger::TDCSector_t::TDCSector1, convert_V1190_Tagger1},
                                        {detector::Tagger::TDCSector_t::TDCSector2, convert_V1190_Tagger2},
                                        {detector::Tagger::TDCSector_t::TDCSector3, convert_V1190_Tagger3}
                                      },
                                      -325, // default offset in ns
                                      std::make_shared<calibration::gui::FitGaus>(),
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
                                               timecuts ? interval<double>{-12, 12} : no_timecut,
                                               std::make_shared<calibration::gui::FitGaus>()
                                               );

    AddCalibration<calibration::CB_Energy>(CB, calibrationDataManager, convert_GeSiCa_SADC,
                                           std::vector<double>{0},    // default pedestal
                                           std::vector<double>{0.07}, // default gain
                                           // default threshold, only used on MC
                                           std::vector<double>{thresholds ? 1.2 : 0.0},
                                           std::vector<double>{1.0}   // default relative gain
                                           );

    AddCalibration<calibration::PID_Energy>(PID, calibrationDataManager, convert_MultiHit16bit,
                                            std::vector<double>{pedestals ? 100.0 : 0.0},   // default pedestals
                                            std::vector<double>{0.014},   // default gain
                                            std::vector<double>{thresholds ? 15.0 : 0.0}, // default Raw threshold
                                            std::vector<double>{0.1},                     // default MC MeV threshold
                                            std::vector<double>{1.0}      // default relative gain
                                            );

    AddCalibration<calibration::TAPS_Energy>(TAPS, calibrationDataManager, convert_MultiHit16bit,
                                             std::vector<double>{pedestals ? 100.0 : 0.0}, // default pedestal
                                             std::vector<double>{0.3}, // default gain
                                             (thresholds ? 5.0 : -std_ext::inf), 0, // default Raw thresholds BaF2/PbWO4
                                             std::vector<double>{thresholds ? 3.4 : 0.0}, // default MC MeV threshold
                                             std::vector<double>{1.0}  // default relative gain
                                             );

    AddCalibration<calibration::TAPS_ShortEnergy>(TAPS, calibrationDataManager, convert_MultiHit16bit, std::vector<double>{pedestals ? 100.0 : 0.0});

    AddCalibration<calibration::TAPSVeto_Energy>(TAPSVeto, calibrationDataManager, convert_MultiHit16bit, pedestals ? 100.0 : 0.0);

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

    //No dedicated smearing for 2017_12 exist yet
    //AddCalibration<calibration::ClusterSmearing>(CB,   "ClusterSmearing",  calibration::ClusterCorrection::Filter_t::MC, calibrationDataManager);
    // ECorr, should be applied after MC smearing, no dedicated ECorr for 2017_12 exists yet
    //AddCalibration<calibration::ClusterECorr>(CB,   "ClusterECorr",  calibration::ClusterCorrection::Filter_t::Both, calibrationDataManager);

    // Prompt random windows - not adapted at all
    AddPromptRange({-20, 20});
    AddRandomRange({ -100,  -50});
    AddRandomRange({  50,   100});
}


bool Setup_2017_12::Matches(const ant::TID& tid) const
{
    if(!std_ext::time_between(tid.Timestamp, "2017-11-28", "2017-12-18"))
        return false;
    return true;
}

double Setup_2017_12::GetElectronBeamEnergy() const {
    return 883.0;
}

void Setup_2017_12::BuildMappings(std::vector<ant::UnpackerAcquConfig::hit_mapping_t>& hit_mappings, std::vector<ant::UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{
    // build the mappings from the given detectors
    // that should provide sane and correct defaults
    Setup::BuildMappings(hit_mappings, scaler_mappings);

    // now you may tweak the mapping at this location here
}

Setup_traits::candidatebuilder_config_t Setup_2017_12::GetCandidateBuilderConfig() const
{
    candidatebuilder_config_t conf;
    conf.PID_Phi_Epsilon = std_ext::degree_to_radian(2.0);
    conf.CB_ClusterThreshold = 12;
    conf.TAPS_ClusterThreshold = 12;
    return conf;
}

Setup_traits::triggersimu_config_t Setup_2017_12::GetTriggerSimuConfig() const
{
    triggersimu_config_t conf;
    conf.Type = triggersimu_config_t::Type_t::CBESum;
    // from https://github.com/padlarson/a2GoAT/blob/AdlarsonAnalysis/src/AdlarsonPhysics.cc#L1018
    // First guesses of edge and width, neither has been extracted from data
    conf.CBESum_Edge = 40; // MeV
    conf.CBESum_Width = 20; // MeV
    return conf;
}

ant::UnpackerA2GeantConfig::promptrandom_config_t Setup_2017_12::GetPromptRandomConfig() const {
    ant::UnpackerA2GeantConfig::promptrandom_config_t conf;
    // default constructed conf has everything disabled
    if(MCTaggerHits) {
        conf.RandomPromptRatio = 0.22; // per unit time interval
        conf.PromptSigma = 0.87;       // in ns
        conf.TimeWindow = {-120, 120};
        conf.PromptOffset = -0.37;
    }
    return conf;
}

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2017_12)


