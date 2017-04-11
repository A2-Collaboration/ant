#include "Setup.h"

#include "Setup_2007_Base.h"

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


struct ImprovedTimeFct2007 : calibration::gui::FitGausPol0 {

    const double window_size; //ns around max

    ImprovedTimeFct2007(const double window = 50.0):
        calibration::gui::FitGausPol0(),
        window_size(window) {}

    ~ImprovedTimeFct2007();

    virtual void SetDefaults(TH1* hist) override {

        if(hist) {
            func->SetParameter(0,hist->GetMaximum());
            double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
            func->SetParameter(1,max_pos);
            SetRange({max_pos - window_size/2.0, max_pos + window_size/2.0});
            func->SetParameter(2,10);

            const auto left  = hist->GetBinContent(hist->FindBin(max_pos-50));
            const auto right = hist->GetBinContent(hist->FindBin(max_pos+50));
            func->SetParameter(3,(left+right)/2.0);
        } else {
            func->SetParameter(0,100);
            func->SetParameter(1,100);
        }

        func->SetParLimits(0,   0, 1E+5);  // positive amplitude
        func->SetParLimits(2, 2.0, 100.0); // sigma

    }
};

ImprovedTimeFct2007::~ImprovedTimeFct2007()
{}

Setup_2007_Base::Setup_2007_Base(const std::string& name, OptionsPtr opt) : Setup(name, opt),
    MCTaggerHits(opt->Get<bool>("MCTaggerHits",false)),
    Trigger(make_shared<detector::Trigger_2007>()),
    Tagger(make_shared<detector::Tagger_2007>()),
    CB(make_shared<detector::CB>()),
    PID(make_shared<detector::PID_2007>()),
    TAPS(make_shared<detector::TAPS_2007>(false, false)), // no cherenkov, don't use sensitive channels
    TAPSVeto(make_shared<detector::TAPSVeto_2007>(false)) // no cherenkov
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
        AddCalibration<calibration::Time>(Tagger,
                                          calibrationDataManager,
                                          convert_CATCH_Tagger,
                                          -325, // default offset in ns
                                          std::make_shared<ImprovedTimeFct2007>(100.0),
                                          timecuts ? interval<double>{-120, 120} : no_timecut
                                          );
        AddCalibration<calibration::Time>(CB,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,      // default offset in ns
                                          std::make_shared<calibration::gui::CBPeakFunction>(),
                                          // before timewalk correction
                                          timecuts ? interval<double>{-20, 200} : no_timecut
                                          );
        AddCalibration<calibration::Time>(PID,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -325,
                                          std::make_shared<calibration::gui::FitGaus>(),
                                          timecuts ? interval<double>{-20, 20} : no_timecut
                                          );
        AddCalibration<calibration::TAPSVeto_Time>(TAPSVeto,
                                                   calibrationDataManager,
                                                   convert_MultiHit16bit,   // for BaF2
                                                   nullptr,                 // for PbWO4
                                                   timecuts ? interval<double>{-12, 12} : no_timecut,
                                                   no_timecut,
                                                   std::make_shared<ImprovedTimeFct2007>(30.0)
                                                   );

        AddCalibration<calibration::CB_Energy>(CB, calibrationDataManager, convert_GeSiCa_SADC,
                                               std::vector<double>{0.0},    // default pedestal
                                               std::vector<double>{0.07},   // default gain
                                               std::vector<double>{thresholds ? 2.0 : 0.0},    // default MeV threshold
                                               std::vector<double>{1.0}   // default relative gain
                                               );

        // CB timing needs timewalk correction
        AddCalibration<calibration::CB_TimeWalk>(CB, calibrationDataManager,
                                                 timecuts ? interval<double>{-10, 20} : no_timecut);

        AddCalibration<calibration::PID_Energy>(PID,
                                                calibrationDataManager,
                                                convert_MultiHit16bit,
                                                std::vector<double>{50,53,30,30,27,45,30,16,40,37,39,33,24,42,38,70,20,23,43,40,19,28,38,37}, /* default pedestals from Acqu */
                                                std::vector<double>{0.014},   // default gain
                                                std::vector<double>{thresholds ? 15.0 : 0.0}, // default Raw threshold
                                                std::vector<double>{0.1},                     // default MC MeV threshold
                                                std::vector<double>{1.0}      // default relative gain
                                                );

        // the PID calibration is a physics module only
        AddCalibration<calibration::PID_PhiAngle>(PID, calibrationDataManager);

        AddCalibration<calibration::TAPSVeto_Energy>(TAPSVeto, calibrationDataManager, convert_MultiHit16bit,
                                                     100.0, 0.010, 0.050, 0.1, 2.5); // increased default gain

        AddCalibration<calibration::TAPS_Time>(TAPS,
                                               calibrationDataManager,
                                               convert_MultiHit16bit,   // for BaF2
                                               nullptr, // for PbWO4
                                               timecuts ? interval<double>{-15, 15} : no_timecut, // for BaF2
                                               no_timecut,
                                               std::make_shared<ImprovedTimeFct2007>(30.0)
                                               );

        AddCalibration<calibration::TAPS_Energy>(TAPS, calibrationDataManager, convert_MultiHit16bit,
                                                 std::vector<double>{100.0}, // default pedestal
                                                 std::vector<double>{0.3}, // default gain
                                                 (thresholds ? 5.0 : -std_ext::inf), 0, // default Raw thresholds BaF2/PbWO4
                                                 std::vector<double>{thresholds ? 1.0 : 0.0}, // default MC MeV threshold
                                                 std::vector<double>{1.0}  // default relative gain
                                                 );

        AddCalibration<calibration::TAPS_ShortEnergy>(TAPS, calibrationDataManager, convert_MultiHit16bit );

    }

double Setup_2007_Base::GetElectronBeamEnergy() const {
        return 1508.0;
}

Setup_traits::candidatebuilder_config_t Setup_2007_Base::GetCandidateBuilderConfig() const {
        candidatebuilder_config_t conf;
        conf.PID_Phi_Epsilon = std_ext::degree_to_radian(2.0);
        conf.CB_ClusterThreshold = 15;
        conf.TAPS_ClusterThreshold = 20;
        return conf;
    }


//AUTO_REGISTER_SETUP(Setup_2007_Base)


}}} // namespace ant::expconfig::setup
