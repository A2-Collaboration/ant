#include "FitVetoBand.h"

#include "BaseFunctions.h"

#include "base/TF1Ext.h"

#include "base/interval.h"
#include "base/Logger.h"

#include "TF1.h"
#include "TH1.h"

#include <algorithm>

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;

void FitVetoBand::Sync()
{
    signal->SetParameters(&(func->GetParameters()[0]));
    bg->SetParameters(    &(func->GetParameters()[3]));
    setRange(signal,GetRange());
    setRange(bg,GetRange());
}

FitVetoBand::FitVetoBand()
{
    signal = functions::exponential::getTF1();
    signal->SetLineColor(kRed);

    signal = new TF1("", "[0]*pow([1], [2]-x)", 0, 1000);

    bg = functions::pol<0>::getTF1();
    bg->SetLineColor(kBlue);

    func = new TF1("","[0]*pow([1], [2]-x)+pol0(3)", 0, 1000);
    func->SetLineColor(kGreen);
    func->SetNpx(1000);

    func->SetParName(0, "Amplitude");
    func->SetParName(1, "base");
    func->SetParName(2, "x_shift");
    func->SetParName(3, "Offset");

    struct AmplitudeKnob : KnobsTF1::TransformedParameterKnob {
        AmplitudeKnob(TF1* Func) :
            TransformedParameterKnob(Func->GetParName(1), Func, 1,
                                     [] (double a, TF1*) { return a*100-99; },
                                     [] (double a, TF1*) { return (a+99)/100; },
                                     IndicatorProperties::Type_t::slider_horizontal
                                     )
        {}
    };
    //AddKnob<AmplitudeKnob>(func);

    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(0), func, 0, IndicatorProperties::Type_t::slider_horizontal);
    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(2), func, 2, IndicatorProperties::Type_t::slider_vertical);
    //AddKnob<KnobsTF1::FixedParameterKnob>(func->GetParName(3), func, 3, IndicatorProperties::Type_t::slider_horizontal);
    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(3), func, 3, IndicatorProperties::Type_t::slider_horizontal);
    AddKnob<KnobsTF1::FixedRangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::FixedRangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
}

FitVetoBand::~FitVetoBand()
{
    delete signal;
    delete bg;
    delete func;
}

void FitVetoBand::Draw()
{
    signal->Draw("same");
    bg->Draw("same");
    func->Draw("same");
}

void FitVetoBand::Fit(TH1 *hist)
{
    EnsureParameterLimits();
    FitFunction::doFit(hist);
    Sync();
}

void gui::FitVetoBand::FitSignal(TH1* hist)
{
    const auto fixedPars = {3};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    Sync();
    UnFixParameters(func, fixedPars);
}

void gui::FitVetoBand::FitBackground(TH1* hist)
{
    const auto fixedPars = {0,1,2};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    Sync();
    UnFixParameters(func, fixedPars);
}

void FitVetoBand::SetDefaults(TH1 *hist)
{
    EnsureParameterLimits();
    func->SetParameter(0, 2);
    func->SetParameter(1, 1.006);
    func->SetParameter(2, 100);
    func->SetParameter(3, hist->GetMinimum(.1));

    Sync();
}

void FitVetoBand::EnsureParameterLimits()
{
    func->SetParLimits(0, 0, 10);
    func->SetParLimits(1, .9, 1.1);
    func->SetParLimits(2, 50, 250);
    func->SetParLimits(3, .5, 5);
}

void FitVetoBand::SetRange(ant::interval<double> i)
{
    setRange(func, i);
    setRange(signal, i);
    setRange(bg, i);
}

ant::interval<double> FitVetoBand::GetRange() const
{
    return getRange(func);
}

FitFunction::SavedState_t FitVetoBand::Save() const
{
    SavedState_t params;
    params.reserve(2+func->GetNpar());
    saveTF1(func,params);

    return params;
}

void FitVetoBand::Load(const SavedState_t &data)
{
    if(data.size() != std::size_t(2+func->GetNpar())) {
        LOG(WARNING) << "Can't load parameters";
        return;
    }

    SavedState_t::const_iterator pos = data.begin();
    loadTF1(pos, func);
    SetRange(getRange(func));

    Sync();

}

double FitVetoBand::Eval(const double energy) const
{
    // make sure energy is not less than miminum of fit range
    double e = energy <= GetRange().Start() ? GetRange().Start() : energy;

    return func->Eval(e);
}

double FitVetoBand::EvalReference(const double energy) const
{
    // make sure energy is not less than miminum of fit range
    double e = energy <= GetRange().Start() ? GetRange().Start() : energy;

    TF1 f(*func);
    // parameters determined from MC, suggested usage range above 100MeV
    f.SetParameters(
                1.968,   // amplitude
                1.0064,  // base
                101.42,  // x shift
                1        // offset
                );

    return f.Eval(e);
}
