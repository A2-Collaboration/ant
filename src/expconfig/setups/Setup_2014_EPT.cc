#include "Setup_2014_EPT.h"

#include "base/std_ext/math.h"
#include "base/Logger.h"

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
#include "calibration/modules/TAPS_ToF.h"
#include "calibration/modules/TAPSVeto_Energy.h"
#include "calibration/modules/TAPSVeto_Time.h"
#include "calibration/modules/Tagger_QDC.h"
#include "calibration/modules/TaggEff.h"
#include "calibration/modules/ClusterSmearing.h"

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


Setup_2014_EPT::Setup_2014_EPT(const string& name, OptionsPtr opt) :
    Setup(name, opt),
    MCTaggerHits(opt->Get<bool>("MCTaggerHits",false))
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
    /// \todo check if 16bit unsigned is correct for all those detectors
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

    // Tagger/EPT QDC measurements need some simple hook
    AddHook<calibration::Tagger_QDC>(EPT->Type, convert_MultiHit16bit);

    // Tagging efficiencies are loaded via a calibration module
    AddCalibration<calibration::TaggEff>(EPT, calibrationDataManager);

    const bool timecuts = !opt->Get<bool>("DisableTimecuts");
    interval<double> no_timecut(-std_ext::inf, std_ext::inf);
    if(!timecuts)
        LOG(INFO) << "Disabling timecuts";

    const bool thresholds = !opt->Get<bool>("DisableThresholds");
    if(!thresholds)
        LOG(INFO) << "Disabling thresholds";

    // then we add the others, and link it to the converters
    AddCalibration<calibration::Time>(EPT,
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
    AddCalibration<calibration::TAPS_Time>(taps,
                                           calibrationDataManager,
                                           convert_MultiHit16bit,   // for BaF2
                                           convert_V1190_TAPSPbWO4, // for PbWO4
                                           timecuts ? interval<double>{-15, 15} : no_timecut, // for BaF2
                                           timecuts ? interval<double>{-25, 25} : no_timecut  // for PbWO4
                                           );
    AddCalibration<calibration::TAPSVeto_Time>(tapsVeto,
                                               calibrationDataManager,
                                               convert_MultiHit16bit,   // for BaF2
                                               convert_V1190_TAPSPbWO4, // for PbWO4
                                               timecuts ? interval<double>{-12, 12} : no_timecut,
                                               timecuts ? interval<double>{-12, 12} : no_timecut
                                               );

    AddCalibration<calibration::CB_Energy>(cb, calibrationDataManager, convert_GeSiCa_SADC,
                                           std::vector<double>{0},    // default pedestal
                                           std::vector<double>{0.07}, // default gain
                                           std::vector<double>{thresholds ? 2.0 : 0.0},    // default threshold
                                           std::vector<double>{1.0}   // default relative gain
                                           );

    AddCalibration<calibration::PID_Energy>(pid, calibrationDataManager, convert_MultiHit16bit );

    AddCalibration<calibration::TAPS_Energy>(taps, calibrationDataManager, convert_MultiHit16bit,
                                             std::vector<double>{100}, // default pedestal
                                             std::vector<double>{0.3}, // default gain
                                             std::vector<double>{thresholds ? 1.0 : 0.0},   // default threshold
                                             std::vector<double>{1.0}  // default relative gain
                                             );

    AddCalibration<calibration::TAPS_ShortEnergy>(taps, calibrationDataManager, convert_MultiHit16bit );

    AddCalibration<calibration::TAPSVeto_Energy>(tapsVeto, calibrationDataManager, convert_MultiHit16bit);

    // enable TAPS shower correction, which is a hook running on list of clusters
    AddCalibration<calibration::TAPS_ShowerCorrection>();

    // add ToF timing to TAPS clusters
    AddCalibration<calibration::TAPS_ToF>(taps, calibrationDataManager);

    // the PID calibration is a physics module only
    AddCalibration<calibration::PID_PhiAngle>(pid, calibrationDataManager);

    // CB timing needs timewalk correction
    AddCalibration<calibration::CB_TimeWalk>(cb, calibrationDataManager,
                                             timecuts ? interval<double>{-10, 20} : no_timecut);

    //Cluster Smearing, Energy. Only activates if root file with histogram present in calibration data folder.
    //Place a file in the MC folder to use MC smearing. Do not put one in the "Data" calibration folder unless
    //you want to smear data as well (probably not...)
    AddCalibration<calibration::ClusterSmearing>(cb,   calibrationDataManager);
    AddCalibration<calibration::ClusterSmearing>(taps, calibrationDataManager);
}

double Setup_2014_EPT::GetElectronBeamEnergy() const {
    return 1604.0;
}

void Setup_2014_EPT::BuildMappings(std::vector<ant::UnpackerAcquConfig::hit_mapping_t>& hit_mappings, std::vector<ant::UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{
    // build the mappings from the given detectors
    // that should provide sane and correct defaults
    Setup::BuildMappings(hit_mappings, scaler_mappings);

    // now you may tweak the mapping at this location here
    // for example, ignore elements
}

ant::ExpConfig::Setup::candidatebuilder_config_t Setup_2014_EPT::GetCandidateBuilderConfig() const
{
    candidatebuilder_config_t conf;
    conf.PID_Phi_Epsilon = std_ext::degree_to_radian(2.0);
    conf.CB_ClusterThreshold = 15;
    conf.TAPS_ClusterThreshold = 20;
    return conf;
}

ant::UnpackerA2GeantConfig::promptrandom_config_t Setup_2014_EPT::GetPromptRandomConfig() const {
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
