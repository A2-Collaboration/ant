#include "TAPSVeto_Time.h"

#include "calibration/fitfunctions/FitGausPol0.h"

#include "expconfig/detectors/TAPSVeto.h"

#include "base/Logger.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

TAPSVeto_Time::TAPSVeto_Time(std::shared_ptr<expconfig::detector::TAPSVeto> tapsveto,
                     shared_ptr<DataManager> calmgr,
                     Calibration::Converter::ptr_t converter_BaF2,
                     Calibration::Converter::ptr_t converter_PbWO4,
                     const interval<double>& timeWindow_BaF2, // default {-inf, inf}
                     const interval<double>& timeWindow_PbWO4, // default {-inf, inf}
                     std::shared_ptr<calibration::gui::PeakingFitFunction> fct
                     ) :
    Time(tapsveto,
         calmgr,
         converter_BaF2, // for BaF2
         -100, // for BaF2
         fct,
         timeWindow_BaF2, // for BaF2
         -0.05 // for BaF2
         )
{
    // DefaultGains, DefaultOffsets and TimeWindows are different for
    // PbWO4 elements

    for(unsigned ch=0; ch<tapsveto->GetNChannels(); ch++)
    {
        if(tapsveto->IsPbWO4(ch)) {
            Converters[ch] = converter_PbWO4;
            DefaultOffsets[ch] = -550;
            DefaultGains[ch] = 1.0;
            TimeWindows[ch] = timeWindow_PbWO4;
        }
    }
}

std::shared_ptr<gui::PeakingFitFunction> TAPSVeto_Time::getDefaultFitFct()
{
    return make_shared<calibration::gui::FitGausPol0>();
}
