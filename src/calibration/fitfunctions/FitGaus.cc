#include "FitGaus.h"

#include "base/interval.h"
#include "base/Logger.h"
#include "BaseFunctions.h"

#include "TF1.h"
#include "TH1.h"

#include <algorithm>

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;

FitGaus::FitGaus()
{
    func = functions::gaus::getTF1();
    func->SetNpx(1000);

    func->SetParName(0,"A");
    func->SetParName(1,"x_{0}");
    func->SetParName(2,"#sigma");

    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(0), func, 0, IndicatorProperties::Type_t::slider_horizontal);
    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(1), func, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::ReferenceParameterKnob>(func->GetParName(2), func, 2, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::RangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
}

FitGaus::~FitGaus()
{
}

void FitGaus::Draw()
{
    func->Draw("same");
}

void FitGaus::Fit(TH1 *hist)
{
    FitFunction::doFit(hist, func);
}

void FitGaus::SetDefaults(TH1 *hist)
{
    if(hist) {
        func->SetParameter(0,hist->GetMaximum());
        const double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
        func->SetParameter(1,max_pos);
        const double sigma = hist->GetRMS();
        func->SetParameter(2, sigma);
        SetRange({max_pos-4*sigma, max_pos+4*sigma});
    } else {
        SetRange({0,200});
        func->SetParameter(0,0.8);
        func->SetParameter(1,100);
        func->SetParameter(2,20);
    }
}

void FitGaus::SetRange(ant::interval<double> i)
{
    setRange(func, i);
}

ant::interval<double> FitGaus::GetRange() const
{
    return getRange(func);
}

FitFunction::SavedState_t FitGaus::Save() const
{
    std::vector<double> params;

    saveTF1(func,params);

    return params;
}

void FitGaus::Load(const SavedState_t &data)
{
    if(data.size() != std::size_t(2+func->GetNpar())) {
        LOG(WARNING) << "Can't load parameters";
        return;
    }
    SavedState_t::const_iterator pos = data.begin();
    loadTF1(pos, func);

    sync();

}

double FitGaus::GetPeakPosition() const
{
    return func->GetParameter(1);
}

double FitGaus::GetPeakWidth() const
{
    return func->GetParameter(2);
}
