#include "Setup_2014_EPT.h"

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


using namespace ant::expconfig;
using namespace ant::expconfig::setup;


Setup_2014_EPT::Setup_2014_EPT(const string& name, SetupOptPtr opt) :
    Setup(name, opt)
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

    bool timecuts = !IsFlagSet("DisableTimecuts");
    interval<double> no_timecut(-std_ext::inf, std_ext::inf);

    AddCalibration<calibration::Scaler>(Detector_t::Type_t::EPT,
                                        convert_ScalerFrequency_Beampolmon);

    // then we add the others, and link it to the converters
    AddCalibration<calibration::Time>(EPT,
                                      calibrationDataManager,
                                      convert_CATCH_Tagger,
                                      -325, // default offset in ns
                                      std::make_shared<calibration::gui::FitGausPol0>(),
                                      timecuts ? interval<double>{-100, 100} : no_timecut
                                      );
    AddCalibration<calibration::Time>(cb,
                                      calibrationDataManager,
                                      convert_CATCH_CB,
                                      -325,      // default offset in ns
                                      std::make_shared<calibration::gui::FitGaus>(),
                                      timecuts ? interval<double>{-10, 100} : no_timecut
                                      );
    AddCalibration<calibration::Time>(pid,
                                      calibrationDataManager,
                                      convert_CATCH_CB,
                                      -325,
                                      std::make_shared<calibration::gui::FitGaus>(),
                                      timecuts ? interval<double>{-7, 7} : no_timecut
                                      );
    AddCalibration<calibration::TAPS_Time>(taps,
                                           calibrationDataManager,
                                           convert_MultiHit16bit,
                                           timecuts ? interval<double>{-10, 10} : no_timecut, // for BaF2
                                           timecuts ? interval<double>{-20, 20} : no_timecut  // for PbWO4
                                           );
    AddCalibration<calibration::Time>(tapsVeto,
                                      calibrationDataManager,
                                      convert_MultiHit16bit,
                                      -100,
                                      std::make_shared<calibration::gui::FitGausPol0>(),
                                      timecuts ? interval<double>{-7, 7} : no_timecut,
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
    AddCalibration<calibration::CB_TimeWalk>(cb, calibrationDataManager,
                                             timecuts ? interval<double>{-5, 10} : no_timecut);

}

double Setup_2014_EPT::GetElectronBeamEnergy() const {
    return 1604.0;
}

bool Setup_2014_EPT::Matches(const ant::THeaderInfo& header) const {
    if(!Setup::Matches(header))
        return false;
    return true;
}

void Setup_2014_EPT::BuildMappings(std::vector<ant::UnpackerAcquConfig::hit_mapping_t>& hit_mappings, std::vector<ant::UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{
    // build the mappings from the given detectors
    // that should provide sane and correct defaults
    Setup::BuildMappings(hit_mappings, scaler_mappings);

    // now you may tweak the mapping at this location here
    // for example, ignore elements
}

ant::ExpConfig::Reconstruct::candidatebuilder_config_t Setup_2014_EPT::GetCandidateBuilderConfig() const
{
    candidatebuilder_config_t conf;
    conf.PID_Phi_Epsilon = std_ext::degree_to_radian(2.0);
    conf.CB_ClusterThreshold = 15;
    conf.TAPS_ClusterThreshold = 20;
    return conf;
}
