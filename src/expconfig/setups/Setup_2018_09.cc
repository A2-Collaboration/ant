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
 * @brief Ant Setup for the Septmeber 2018 beamtime using the new Tagger
 */
class Setup_2018_09 : public Setup_2017Plus_NewTagger_Base
{
protected:
    const std::shared_ptr<detector::Tagger_2018_03> Tagger;

public:

    Setup_2018_09(const string& name, OptionsPtr opt) :
        Setup_2017Plus_NewTagger_Base(name, opt),
        Tagger(make_shared<detector::Tagger_2018_03>())
    {
        // add the specific Tagger cabling
        AddDetector(Tagger);



        // Broken, BadTDC or NoCalib elements
        CB->SetElementFlag(Detector_t::ElementFlag_t::BadTDC, {17,265,582,586,672,678,696});
        CB->SetElementFlag(Detector_t::ElementFlag_t::Broken, {99,125,162,461,503});
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::Broken, TAPS->GetPbWO4Channels()); //All the PbWO were turned off
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::Broken, {114,137}); //And two more
        TAPS->SetElementFlag(Detector_t::ElementFlag_t::NoCalibFill, {138,145,218,283,346,356,357,364});
        TAPSVeto->SetElementFlag(Detector_t::ElementFlag_t::Broken, {0,1,2,3,5,6,15,31,36,41,64,65,66,71,75,92,119,128,129,
                                                                     130,188,191,192,193,194,195,203,242,243,253,254,256,257,
                                                                     258,263,287,288,291,292,307,311,320,321,322,337,349,350,354});
        Tagger->SwitchOffElementRange(0, 47);

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
                                                                                                                     Trigger->Reference_V1190_TaggerTDC3_2,
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
                                                        interval<double>{-300, 300} :
                                                        interval<double>{-std_ext::inf, std_ext::inf}
                                                        );

    }


    // change beam energy (if not 883 MeV) like this
    //double GetElectronBeamEnergy() const override {
    //    return 883.0;
    //}

    bool Matches(const ant::TID& tid) const override
    {
        if(!std_ext::time_between(tid.Timestamp, "2018-09-11", "2018-10-01"))
            return false;
        return true;
    }

    triggersimu_config_t GetTriggerSimuConfig() const override
    {
        auto conf = Setup_2017Plus_NewTagger_Base::GetTriggerSimuConfig();
        return conf;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2018_09)


}}} // namespace ant::expconfig::setup
