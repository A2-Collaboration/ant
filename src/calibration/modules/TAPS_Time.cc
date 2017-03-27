#include "TAPS_Time.h"

#include "calibration/fitfunctions/FitGausPol0.h"

#include "expconfig/detectors/TAPS.h"

#include "base/Logger.h"

#include <cstdint>

#include "TH1.h"
#include "TF1.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;


class TAPSTimeFunction : public calibration::gui::FitGausPol0 {
public:
    using FitGausPol0::FitGausPol0;

    virtual void SetDefaults(TH1* hist) override;

};


TAPS_Time::TAPS_Time(shared_ptr<expconfig::detector::TAPS> taps,
                     shared_ptr<DataManager> calmgr,
                     Calibration::Converter::ptr_t converter_BaF2,
                     Calibration::Converter::ptr_t converter_PbWO4,
                     const interval<double>& timeWindow_BaF2, // default {-inf, inf}
                     const interval<double>& timeWindow_PbWO4, // default {-inf, inf}
                     std::shared_ptr<calibration::gui::PeakingFitFunction> fct
                     ) :
    Time(taps,
         calmgr,
         converter_BaF2, // for BaF2
         -170, // for BaF2
         fct,
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


void TAPSTimeFunction::SetDefaults(TH1 *hist)
{

    if(hist) {

        const auto height = hist->GetMaximum();

        const double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
        SetRange({max_pos - 20, max_pos + 20});
        const auto range = GetRange();
        const auto height_right = hist->GetBinContent(hist->FindBin(range.Start()));
        const auto height_left = hist->GetBinContent(hist->FindBin(range.Stop()));
        const auto offset = (height_left+height_right)/2.0;

        func->SetParameter(0, height - offset);
        func->SetParLimits(0, 0.0, 2 * height); // positive amplitude

        func->SetParameter(1,max_pos);
        func->SetParLimits(1,range.Start(), range.Stop());

        func->SetParameter(2,10.0);
        func->SetParLimits(2, 1.0, GetRange().Length() / 2.0);

        func->SetParameter(3, offset);
        func->SetParLimits(3, 0.0, height);
    } else {
        func->SetParameter(0,100.0);
        func->SetParameter(1,100.0);
        func->SetParameter(2,2.0);
        func->SetParameter(3,0.0);
    }
}

std::shared_ptr<gui::PeakingFitFunction> TAPS_Time::getDefaultFitFct()
{
    return make_shared<TAPSTimeFunction>();
}


