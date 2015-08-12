#include "FitGaus.h"

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

FitGaus::FitGaus()
{
    func = functions::gaus::getFT1();
    func->SetNpx(1000);
    Addknob<KnobsTF1::ParameterKnob>("A",     func, 0, GUIElementDescription::GUI_Type::slider_horizontal, kBlue, 3);
    Addknob<KnobsTF1::ParameterKnob>("x_{0}", func, 1, GUIElementDescription::GUI_Type::slider_vertical,   kBlue, 3);
    Addknob<MyWKnob>("#sigma",func);
    Addknob<KnobsTF1::RangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    Addknob<KnobsTF1::RangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
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
    hist->Fit(func,"RBQN");
}

void FitGaus::SetDefaults(TH1 *hist)
{
    SetRange({0,400});
    if(hist) {
        func->SetParameter(0,hist->GetMaximum());
    } else {
        func->SetParameter(0,1000);
    }
    func->SetParameter(1,135);
    func->SetParameter(2,20);
}

void FitGaus::SetRange(ant::interval<double> i)
{
    setRange(func, i);
}

ant::interval<double> FitGaus::GetRange() const
{
    return getRange(func);
}

void FitGaus::SetPoints(int n)
{
    func->SetNpx(n);
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
        LOG(WARNING) << "Can't load parametes";
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

FitGaus::MyWKnob::MyWKnob(const std::string &n, TF1 *Func):
    VirtualKnob(n,{GUIElementDescription::GUI_Type::slider_vertical,kBlue,3}),
    func(Func)
{
}

double FitGaus::MyWKnob::get() const
{
    return func->GetParameter(1) + func->GetParameter(2);
}

void FitGaus::MyWKnob::set(double a)
{
    auto v = a - func->GetParameter(1);
    func->SetParameter(2,v);
}
