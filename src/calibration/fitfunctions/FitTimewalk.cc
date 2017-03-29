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
    FitFunction::doFit(hist);
}

void FitTimewalk::FitSignal(TH1* hist)
{
    using p = functions::timewalk::p;
    const auto fixedPars = {p::Offset,p::Slope,p::E0}; // keep offset, slope, edge=E_0
    EnsureParameterLimits();
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    UnFixParameters(func, fixedPars);
}

void FitTimewalk::FitBackground(TH1* hist)
{
    using p = functions::timewalk::p;
    const auto fixedPars = {p::Scale,p::Pow,p::Exp}; // keep scale, power, exponent
    EnsureParameterLimits();
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    UnFixParameters(func, fixedPars);
}

void FitTimewalk::SetDefaults(TH1*)
{
    using p = functions::timewalk::p;
    func->SetParameter(p::Offset,    0); // Offset
    func->SetParameter(p::Scale,    50); // scale
    func->SetParameter(p::E0,      1.1); // E_0
    func->SetParameter(p::Pow,     1); // power
    func->SetParameter(p::Exp,     0.5); // exp scale
    func->SetParameter(p::Slope,  0); // linear slope
}

void FitTimewalk::EnsureParameterLimits()
{
    using p = functions::timewalk::p;
    func->SetParLimits(p::Offset, -100, 100);
    func->SetParLimits(p::Scale, 0, 1000);
    func->SetParLimits(p::E0, 0.5, 1.3);
    func->SetParLimits(p::Pow, 0.0001, 3);
    func->SetParLimits(p::Exp, 0, 5);
    func->SetParLimits(p::Slope, -0.5, 0);
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
}

double FitTimewalk::Eval(double log10_raw_energy) const
{
    // make sure raw energy is not less that miminum of fit range
    if(log10_raw_energy <= GetRange().Start())
        log10_raw_energy = GetRange().Start();
    // make sure that raw energy is not beyond asymptote defined by E0
    const auto E0 = func->GetParameter(functions::timewalk::p::E0);
    if(log10_raw_energy <= E0)
        // as raw energy is digitized, +1 is meaningful here on non-log scale
        // we get the value "very close" to the asymptote
        log10_raw_energy = std::log10(std::pow(10, E0) + 1);
    return func->Eval(log10_raw_energy);
}


