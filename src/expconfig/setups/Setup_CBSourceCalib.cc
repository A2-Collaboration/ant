#include "Setup.h"

#include "base/Logger.h"

#include "detectors/CB.h"
#include "detectors/Trigger.h"


#include "calibration/modules/Time.h"
#include "calibration/modules/CB_Energy.h"
#include "calibration/modules/CB_SourceCalib.h"
#include "calibration/converters/MultiHit.h"
#include "calibration/converters/CATCH_TDC.h"
#include "calibration/converters/GeSiCa_SADC.h"


using namespace std;

namespace ant {
namespace expconfig {
namespace setup {


class Setup_CBSourceCalib : public Setup
{
public:

    Setup_CBSourceCalib(const std::string& name, OptionsPtr opt) : Setup(name, opt)
    {

        auto cb = make_shared<detector::CB>();
        AddDetector(cb);

        const auto& convert_CATCH_CB = make_shared<calibration::converter::CATCH_TDC>(
                    detector::Trigger::Reference_CATCH_CBCrate
                    );
        const auto& convert_GeSiCa_SADC = make_shared<calibration::converter::GeSiCa_SADC>();
        AddHook(convert_CATCH_CB);

        const bool timecuts = !opt->Get<bool>("DisableTimecuts");
        interval<double> no_timecut(-std_ext::inf, std_ext::inf);
        if(!timecuts)
            LOG(INFO) << "Disabling timecuts";

//        const bool thresholds = !Options->Get<bool>("DisableThresholds");
//        if(!thresholds)
//            LOG(INFO) << "Disabling thresholds";


        AddCalibration<calibration::Time>(cb,
                                          calibrationDataManager,
                                          convert_CATCH_CB,
                                          -500,      // default offset in ns
                                          std::make_shared<calibration::gui::CBPeakFunction>(),
                                          // before timewalk correction
                                          timecuts ? interval<double>{-1000, 1000} : no_timecut
                                          );

        AddCalibration<calibration::CB_SourceCalib>(cb, calibrationDataManager, convert_GeSiCa_SADC);
//                                               0,    // default pedestal
//                                               1, // default gain
//                                               0,     // default threshold
//                                               1.0   // default relative gain
//                                               );
    }

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const override {
        Setup::BuildMappings(hit_mappings, scaler_mappings);
        auto& cb_reftiming = detector::Trigger::Reference_CATCH_CBCrate;
        hit_mappings.emplace_back(
                    cb_reftiming.LogicalChannel,
                    cb_reftiming.AcquRawChannel
                    );
    }
};

AUTO_REGISTER_SETUP(Setup_CBSourceCalib)
}
}
}
