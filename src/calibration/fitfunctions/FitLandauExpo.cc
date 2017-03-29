#include "FitLandauExpo.h"

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

void FitLandauExpo::Sync()
{
    signal->SetParameters(&(func->GetParameters()[0]));
    bg->SetParameters(    &(func->GetParameters()[3]));
    setRange(signal,GetRange());
    setRange(bg,GetRange());
}

FitLandauExpo::FitLandauExpo()
{
    signal = new TF1("","landau", 0, 200);
    signal->SetLineColor(kRed);

    bg = functions::exponential::getTF1();
    bg->SetLineColor(kBlue);

    func = new TF1("","landau+expo(3)", 0, 200);
    func->SetLineColor(kGreen);
    func->SetNpx(1000);

    func->SetParName(0, "A");
    func->SetParName(1, "MPV"); // most probable value
    func->SetParName(2, "#sigma");
    func->SetParName(3, "p_{0}");
    func->SetParName(4, "p_{1}");

    struct AmplitudeKnob : KnobsTF1::TransformedParameterKnob {
        AmplitudeKnob(TF1* Func) :
            TransformedParameterKnob(Func->GetParName(0), Func, 0,
                                     [] (double a, TF1*) { return a/5.5; },
                                     [] (double a, TF1*) { return a*5.5; },
                                     IndicatorProperties::Type_t::slider_horizontal
                                     )
        {}
    };
    AddKnob<AmplitudeKnob>(func);

    // this ensure the 100% correct position of the most probable value,
    // see https://root.cern.ch/root/html524/TMath.html#TMath:Landau
    MPV_trafo = [] (double a, TF1* f) { return a - f->GetParameter(2)*0.22278; };
    MPV_trafo_inverse = [] (double a, TF1* f) { return a + f->GetParameter(2)*0.22278; };
    struct MPVKnob : KnobsTF1::TransformedParameterKnob {
        MPVKnob(TF1* Func, transformation_t trafo, transformation_t trafo_inverse) :
            TransformedParameterKnob(Func->GetParName(1), Func, 1,
                                     trafo,
                                     trafo_inverse,
                                     IndicatorProperties::Type_t::slider_vertical
                                     )
        {}
    };
    AddKnob<MPVKnob>(func, MPV_trafo, MPV_trafo_inverse);

    AddKnob<KnobsTF1::ReferenceParameterKnob>(func->GetParName(2), func, 2, 1, IndicatorProperties::Type_t::slider_vertical);
    //AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(3), func, 3, IndicatorProperties::Type_t::slider_horizontal);
    //AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(4), func, 4, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::RangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
}

FitLandauExpo::~FitLandauExpo()
{
    delete signal;
    delete bg;
    delete func;
}

void FitLandauExpo::Draw()
{
    signal->Draw("same");
    bg->Draw("same");
    func->Draw("same");
}

void FitLandauExpo::Fit(TH1 *hist)
{
    FitFunction::doFit(hist);
    Sync();
}

void gui::FitLandauExpo::FitBackground(TH1* hist)
{
    const auto fixedPars = {0,1,2};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    Sync();
    UnFixParameters(func, fixedPars);
}

void gui::FitLandauExpo::FitSignal(TH1* hist)
{
    const auto fixedPars = {3,4};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    Sync();
    UnFixParameters(func, fixedPars);
}

void FitLandauExpo::SetDefaults(TH1 *hist)
{
    func->SetParameter(2, 3.0);
    SetRange({0, 200});

    if(hist) {
        auto range = GetRange();
        func->SetParameter(0, hist->GetMaximum()/3.0);
        const double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
        func->SetParameter(1,max_pos);
        func->SetParameter(3, log(range.Clip(max_pos)));
    } else {
        func->SetParameter(0,1000);
        func->SetParameter(1,100);
    }
    func->SetParameter(4, -0.7);

    Sync();
}

void FitLandauExpo::SetRange(ant::interval<double> i)
{
    setRange(func, i);
    setRange(signal, i);
    setRange(bg, i);
    // x_0 peak position must be in range
    func->SetParLimits(1, i.Start(), i.Stop());
}

ant::interval<double> FitLandauExpo::GetRange() const
{
    return getRange(func);
}

FitFunction::SavedState_t FitLandauExpo::Save() const
{
    SavedState_t params;
    params.reserve(2+func->GetNpar());
    saveTF1(func,params);

    return params;
}

void FitLandauExpo::Load(const SavedState_t &data)
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

double FitLandauExpo::GetPeakPosition() const
{
    return MPV_trafo(func->GetParameter(1), func);
}

double FitLandauExpo::GetPeakWidth() const
{
    return func->GetParameter(2);
}

double FitLandauExpo::SignalToBackground(const double x) const
{
    const auto s = func->Eval(x);
    const auto b = bg->Eval(x);

    return (s-b)/(s+b);
}
