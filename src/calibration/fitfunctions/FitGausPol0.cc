#include "FitGausPol0.h"

#include "base/interval.h"
#include "base/Logger.h"
#include "TF1Knobs.h"
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
    AddKnob<AmpKnop>("A",func);
    AddKnob<KnobsTF1::ParameterKnob>("x_{0}", func, 1, IndicatorProperties::Type_t::slider_vertical,   kBlue, 3);
    AddKnob<SigmaKnob>("#sigma",func);
    AddKnob<KnobsTF1::ParameterKnob>("c",     func, 3, IndicatorProperties::Type_t::slider_horizontal, kRed, 3);
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
    hist->Fit(func,"RBQN");
}

void FitGausPol0::SetDefaults(TH1 *hist)
{
    if(hist) {
        func->SetParameter(0,hist->GetMaximum());
        double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
        func->SetParameter(1,max_pos);
        SetRange({max_pos-20, max_pos+20});
        func->SetParameter(2,10);
        func->SetParameter(3,(hist->GetMaximum() - hist->GetMinimum()) / 2 );
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

FitGausPol0::SigmaKnob::SigmaKnob(const std::string &n, TF1 *Func):
    IndicatorKnob(n,{IndicatorProperties::Type_t::slider_vertical,kBlue,3}),
    func(Func)
{
}

double FitGausPol0::SigmaKnob::get() const
{
    return func->GetParameter(1) + func->GetParameter(2);
}

void FitGausPol0::SigmaKnob::set(double a)
{
    auto v = a - func->GetParameter(1);
    func->SetParameter(2,v);
}


FitGausPol0::AmpKnop::AmpKnop(const std::string& name, TF1* Func):
    IndicatorKnob(name,{IndicatorProperties::Type_t::slider_horizontal,kBlue,3}),
    func(Func)
{

}

double FitGausPol0::AmpKnop::get() const
{
    return func->GetParameter(0) + func->GetParameter(3);
}

void FitGausPol0::AmpKnop::set(double a)
{
    auto v = a - func->GetParameter(3);
    func->SetParameter(0,v);
}
