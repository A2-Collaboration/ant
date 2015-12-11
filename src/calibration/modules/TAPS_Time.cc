#include "TAPS_Time.h"

#include "calibration/fitfunctions/FitGausPol0.h"

#include "expconfig/detectors/TAPS.h"

#include "base/Logger.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

TAPS_Time::TAPS_Time(shared_ptr<expconfig::detector::TAPS> taps,
                     shared_ptr<DataManager> calmgr,
                     Calibration::Converter::ptr_t converter_BaF2,
                     Calibration::Converter::ptr_t converter_PbWO4,
                     const interval<double>& timeWindow_BaF2, // default {-inf, inf}
                     const interval<double>& timeWindow_PbWO4 // default {-inf, inf}
                     ) :
    Time(taps,
         calmgr,
         converter_BaF2, // for BaF2
         -170, // for BaF2
         std::make_shared<calibration::gui::FitGausPol0>(),
         timeWindow_BaF2, // for BaF2
         -0.100 // for BaF2
         )
{
    // DefaultGains, DefaultOffsets and TimeWindows are different for
    // PbWO4 elements

    for(unsigned ch=0; ch<taps->GetNChannels(); ch++)
    {
        if(taps->IsPbWO4(ch)) {
            Converters[ch] = converter_PbWO4;
            DefaultOffsets[ch] = -550;
            DefaultGains[ch] = 1.0;
            TimeWindows[ch] = timeWindow_PbWO4;
        }
    }
}
