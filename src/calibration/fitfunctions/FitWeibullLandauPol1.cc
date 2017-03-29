#include "FitWeibullLandauPol1.h"

#include "BaseFunctions.h"

#include "base/TF1Ext.h"

#include "base/std_ext/math.h"

#include "base/interval.h"
#include "base/Logger.h"

#include "TF1.h"
#include "TH1.h"

#include <algorithm>

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;

void FitWeibullLandauPol1::Sync()
{
    signal->SetParameters(&(func->GetParameters()[0]));
    bg->SetParameters(    &(func->GetParameters()[5]));
    setRange(signal,GetRange());
    setRange(bg,GetRange());
}

FitWeibullLandauPol1::FitWeibullLandauPol1()
{
    signal = functions::WeibullLandau::getTF1();
    signal->SetLineColor(kRed);

    bg = functions::pol<1>::getTF1();
    bg->SetLineColor(kBlue);

    func = functions::WeibullLandauPol<1>::getTF1();
    func->SetLineColor(kGreen);
    func->SetNpx(1000);

    func->SetParName(0, "#lambda");  // weibull inverse scale parameter
    func->SetParName(1, "k");        // weibull shape parameter
    func->SetParLimits(1, 0.1, 10);
    func->SetParName(2, "A");        // amplitude
    func->SetParName(3, "MPV");      // most probable value
    func->SetParName(4, "#sigma");   // landau sigma
    func->SetParLimits(4, 0.1, 5);
    func->SetParName(5, "b");        // offset
    func->SetParName(6, "m");        // pol1 slope


    // this ensure the 100% correct position of the most probable value,
    // see https://root.cern.ch/root/html524/TMath.html#TMath:Landau
/*    MPV_trafo = [] (double a, TF1* f) { return a - f->GetParameter(4)*0.22278; };
    MPV_trafo_inverse = [] (double a, TF1* f) { return a + f->GetParameter(4)*0.22278; };
    struct MPVKnob : KnobsTF1::TransformedParameterKnob {
        MPVKnob(TF1* Func, transformation_t trafo, transformation_t trafo_inverse) :
            TransformedParameterKnob(Func->GetParName(3), Func, 3,
                                     trafo,
                                     trafo_inverse,
                                     IndicatorProperties::Type_t::slider_vertical
                                     )
        {}
    };
    AddKnob<MPVKnob>(func, MPV_trafo, MPV_trafo_inverse);

    AddKnob<KnobsTF1::ReferenceParameterKnob>(func->GetParName(1), func, 1, IndicatorProperties::Type_t::slider_vertical);  // weibull k
    //AddKnob<KnobsTF1::ReferenceParameterKnob>(func->GetParName(4), func, 4, 1, IndicatorProperties::Type_t::slider_vertical);  // landau sigma
    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(5), func, 5, IndicatorProperties::Type_t::slider_horizontal);  // pol1 y axis
*/
    AddKnob<KnobsTF1::RangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
}

FitWeibullLandauPol1::~FitWeibullLandauPol1()
{
    delete signal;
    delete bg;
    delete func;
}

void FitWeibullLandauPol1::Draw()
{
    signal->Draw("same");
    bg->Draw("same");
    func->Draw("same");
}

void FitWeibullLandauPol1::Fit(TH1 *hist)
{
    FitFunction::doFit(hist);
    Sync();
}

void gui::FitWeibullLandauPol1::FitBackground(TH1* hist)
{
    const auto fixedPars = {0,1,2,3,4};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    Sync();
    UnFixParameters(func, fixedPars);
}

void gui::FitWeibullLandauPol1::FitSignal(TH1* hist)
{
    const auto fixedPars = {5,6};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    Sync();
    UnFixParameters(func, fixedPars);
}

void FitWeibullLandauPol1::SetDefaults(TH1 *hist)
{
    func->SetParameter(1, 0.2);
    func->SetParameter(3, 5.5);
    func->SetParameter(4, 1.1);

    if (hist) {
        const auto range = GetRange();
        func->SetParameter(0, hist->GetMaximum()*3);
        func->SetParameter(2, hist->GetEntries()*hist->GetMaximum());
        func->SetParameter(5, hist->GetMaximum()/10);
        const double y1 = hist->Interpolate(range.Start());
        const double y2 = hist->Interpolate(range.Stop());
        const double slope = (y2 - y1)/(range.Stop() - range.Start());
        func->SetParameter(6, slope);
    } else {
        func->SetParameter(0, 100);
        func->SetParameter(2, 1e6);
        func->SetParameter(6, -20);
    }

    const auto range = GetRange();
    // Landau peak position must be in range
    func->SetParLimits(3, range.Start(), range.Stop());
    // Weibull distribution only defined for positive parameters
    func->SetParLimits(0, 0, 1e10);
    func->SetParLimits(1, 0.1, 10);
    // amplitude of Landau has to be positive
    func->SetParLimits(2, 10, 1e20);
    // width of Landau
    func->SetParLimits(4, 0.1, 5);

    Sync();
}

void FitWeibullLandauPol1::SetRange(ant::interval<double> i)
{
    setRange(func, i);
    setRange(signal, i);
    setRange(bg, i);

}

ant::interval<double> FitWeibullLandauPol1::GetRange() const
{
    return getRange(func);
}

FitFunction::SavedState_t FitWeibullLandauPol1::Save() const
{
    SavedState_t params;
    params.reserve(2+func->GetNpar());
    saveTF1(func,params);

    return params;
}

void FitWeibullLandauPol1::Load(const SavedState_t &data)
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

double FitWeibullLandauPol1::GetWeibullPeak() const
{
    const double lambda = func->GetParameter(0);
    const double k = func->GetParameter(1);
    const double mean = 1/lambda * tgamma(1 + 1/k);
    return mean;
}

double FitWeibullLandauPol1::GetWeibullWidth() const
{
    const double lambda = func->GetParameter(0);
    const double k = func->GetParameter(1);
    const double var = 1/lambda/lambda * ( tgamma(1 + 2/k) - std_ext::sqr(tgamma(1 + 1/k)) );
    return sqrt(var);
}

double FitWeibullLandauPol1::GetLandauPeak() const
{
    return MPV_trafo(func->GetParameter(3), func);
}

double FitWeibullLandauPol1::GetLandauWidth() const
{
    return func->GetParameter(4);
}

double FitWeibullLandauPol1::GetPeakPosition() const
{
    return func->GetMaximumX(2,5);
}

double FitWeibullLandauPol1::GetPeakWidth() const
{
    // calculate FWHM within the defined range for the signal function
    const auto range = GetRange();
    const double low = range.Start(), high = range.Stop();
    const double max = signal->GetMaximum(low, high);
    const double max_pos = signal->GetMaximumX(low, high);
    const double x_low = signal->GetX(max/2, low, max_pos);
    const double x_high = signal->GetX(max/2, max_pos, high);
    return x_high - x_low;
}

double FitWeibullLandauPol1::SignalToBackground(const double x) const
{
    const auto s = func->Eval(x);
    const auto b = bg->Eval(x);

    return (s-b)/(s+b);
}
