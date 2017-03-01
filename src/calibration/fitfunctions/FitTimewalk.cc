#include "FitTimewalk.h"

#include "BaseFunctions.h"

#include "base/interval.h"
#include "base/Logger.h"
#include "base/TF1Ext.h"

#include "TF1.h"
#include "TH1.h"

#include <algorithm>

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;



FitTimewalk::FitTimewalk()
{
    // p[0] + p[5]*x0 +  p[1]*std::exp(-p[4]*x0 - p[3]*std::log(x0));
    func = functions::timewalk::getTF1();
    func->SetNpx(1000);
    func->SetParName(0, "Offset");
    func->SetParName(2, "E_{0}");

    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(0), func, 0, IndicatorProperties::Type_t::slider_horizontal);

    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(2), func, 2, IndicatorProperties::Type_t::slider_vertical);

    AddKnob<KnobsTF1::RangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
}

FitTimewalk::~FitTimewalk()
{
}

void FitTimewalk::Draw()
{
    func->Draw("same");
}

void FitTimewalk::Fit(TH1 *hist)
{
    EnsureParameterLimits();
    FitFunction::doFit(hist, func);
}

void FitTimewalk::FitSignal(TH1* hist)
{
    const auto fixedPars = {0,2,5}; // keep offset, slope, edge=E_0
    EnsureParameterLimits();
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist, func);
    UnFixParameters(func, fixedPars);
}

void FitTimewalk::FitBackground(TH1* hist)
{
    const auto fixedPars = {1,3,4}; // keep scale, power, exponent
    EnsureParameterLimits();
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist, func);
    UnFixParameters(func, fixedPars);
}

void FitTimewalk::SetDefaults(TH1* h)
{
    if(loaded)
        return;

    func->SetParameter(0, 0); // Offset
    func->SetParameter(1, 50);  // scale
    func->SetParameter(2, 20);  // E_0
    func->SetParameter(3, 0.1); // power
    func->SetParameter(4, 0.05); // exp scale
    func->SetParameter(5, -0.01); // linear slope

    SetRange({25, h->GetXaxis()->GetXmax()});

    loaded = true;
}

void FitTimewalk::EnsureParameterLimits()
{
    func->SetParLimits(0, -100, 100);
    func->SetParLimits(1, 0, 1000);
    func->SetParLimits(2, -10, 30);
    func->SetParLimits(3, 0.0001, 3);
    func->SetParLimits(4, 0, 5);
    func->SetParLimits(5, -0.5, 0);
}

void FitTimewalk::SetRange(ant::interval<double> i)
{
    setRange(func, i);
}

ant::interval<double> FitTimewalk::GetRange() const
{
    return getRange(func);
}

FitFunction::SavedState_t FitTimewalk::Save() const
{
    std::vector<double> params;

    saveTF1(func,params);

    return params;
}

void FitTimewalk::Load(const SavedState_t &data)
{

    if(data.size() != std::size_t(2+func->GetNpar())) {
        LOG(WARNING) << "Can't load parameters";
        return;
    }
    SavedState_t::const_iterator pos = data.begin();
    loadTF1(pos, func);
    Sync();
    loaded = true;
}

double FitTimewalk::Eval(double energy)
{
    return func->Eval(energy);
}


