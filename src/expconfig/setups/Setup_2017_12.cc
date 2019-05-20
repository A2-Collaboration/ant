#include "Setup_2017Plus_NewTagger_Base.h"

#include "calibration/modules/Tagger_QDC.h"
#include "calibration/modules/TaggEff.h"
#include "calibration/modules/NewTagger_Time.h"

#include "calibration/converters/MultiHit.h"
#include "calibration/converters/MultiHitReference.h"


using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Ant Setup for the December 2017 beamtime using the new Tagger
 */
class Setup_2017_12 : public Setup_2017Plus_NewTagger_Base
{
protected:
    const std::shared_ptr<detector::Tagger_2017_12> Tagger;

public:

    Setup_2017_12(const string& name, OptionsPtr opt) :
        Setup_2017Plus_NewTagger_Base(name, opt),
        Tagger(make_shared<detector::Tagger_2017_12>())
    {
        SetTimeRange("2017-11-28", "2017-12-18");

        // add the specific Tagger cabling
        AddDetector(Tagger);


        // Broken, BadTDC or NoCalib elements
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {125,461,557,678});
        CB->SetElementFlag(Detector_t::ElementFlag_t::BadTDC, {17,265,461,582,586,672,678,696});
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
                                                                                                                     Trigger->Reference_V1190_TaggerTDC3_1,
                                                                                                                     calibration::converter::Gains::V1190_TDC
                                                                                                                     );


        // the order of the reconstruct hooks is important
        // add both CATCH converters and the V1190 first,
        // since they need to scan the detector read for their reference hit
        AddHook(convert_V1190_Tagger1);
        AddHook(convert_V1190_Tagger2);
        AddHook(convert_V1190_Tagger3);

        // Tagger/EPT QDC measurements need some simple hook
        AddHook<calibration::Tagger_QDC>(Tagger->Type, convert_MultiHit16bit);

        // Tagging efficiencies are loaded via a calibration module
        AddCalibration<calibration::TaggEff>(Tagger, calibrationDataManager);

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
                                                    !opt->Get<bool>("DisableTimecuts") ?
                                                        interval<double>{-120, 120} :
                                                        interval<double>{-std_ext::inf, std_ext::inf}
                                                        );
    }


    triggersimu_config_t GetTriggerSimuConfig() const override
    {
        auto conf = Setup_2017Plus_NewTagger_Base::GetTriggerSimuConfig();
        return conf;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2017_12)


}}} // namespace ant::expconfig::setup
