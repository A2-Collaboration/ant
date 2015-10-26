#include "FitGausPol0.h"

#include "base/interval.h"
#include "base/Logger.h"
#include "BaseFunctions.h"

#include "TF1.h"
#include "TH1.h"

#include <algorithm>

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;

FitGausPol0::FitGausPol0()
{
    func = functions::GausPol<0>::getTF1();
    func->SetNpx(1000);

    func->SetParName(0,"A");
    func->SetParName(1,"x_{0}");
    func->SetParName(2,"#sigma");
    func->SetParName(3,"offset");

    AddKnob<KnobsTF1::ReferenceParameterKnob>(func->GetParName(0), func, 0, 3, IndicatorProperties::Type_t::slider_horizontal);
    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(1), func, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::ReferenceParameterKnob>(func->GetParName(2), func, 2, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(3), func, 3, IndicatorProperties::Type_t::slider_horizontal, kRed, 3);
    AddKnob<KnobsTF1::RangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
}

FitGausPol0::~FitGausPol0()
{
}

void FitGausPol0::Draw()
{
    func->Draw("same");
}

void FitGausPol0::Fit(TH1 *hist)
{
    FitFunction::doFit(hist, func);
}

void FitGausPol0::SetDefaults(TH1 *hist)
{
    if(hist) {
        func->SetParameter(0,hist->GetMaximum());
        double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
        func->SetParameter(1,max_pos);
        SetRange({max_pos-20, max_pos+20});
        func->SetParameter(2,10);
        func->SetParameter(3,0);
    } else {
        func->SetParameter(0,100);
        func->SetParameter(1,100);
    }
}

void FitGausPol0::SetRange(ant::interval<double> i)
{
    setRange(func, i);
}

ant::interval<double> FitGausPol0::GetRange() const
{
    return getRange(func);
}

FitFunction::SavedState_t FitGausPol0::Save() const
{
    std::vector<double> params;

    saveTF1(func,params);

    return params;
}

void FitGausPol0::Load(const SavedState_t &data)
{
    if(data.size() != std::size_t(2+func->GetNpar())) {
        LOG(WARNING) << "Can't load parameters";
        return;
    }
    SavedState_t::const_iterator pos = data.begin();
    loadTF1(pos, func);

    sync();

}

double FitGausPol0::GetPeakPosition() const
{
    return func->GetParameter(1);
}

double FitGausPol0::GetPeakWidth() const
{
    return func->GetParameter(2);
}
