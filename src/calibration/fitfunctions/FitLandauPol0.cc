#include "FitLandauPol0.h"

#include "base/interval.h"
#include "base/Logger.h"
#include "BaseFunctions.h"

#include "TF1.h"
#include "TH1.h"

#include <algorithm>

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;

FitLandauPol0::FitLandauPol0()
{
    func = new TF1("","landau+pol0(3)", 0, 200);
    func->SetNpx(1000);

    func->SetParName(0,"A"); // amplitude
    func->SetParName(1,"MPV"); // most probable value
    func->SetParName(2,"#sigma");
    func->SetParName(3,"Offset");


    struct AmplitudeKnob : KnobsTF1::TransformedParameterKnob {
        AmplitudeKnob(TF1* Func) :
            TransformedParameterKnob(Func->GetParName(0), Func, 0,
                                     [] (double a, TF1* f) { return a/5.5 + f->GetParameter(3); },
                                     [] (double a, TF1* f) { return (a - f->GetParameter(3))*5.5; },
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
    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(3), func, 3, IndicatorProperties::Type_t::slider_horizontal);

    AddKnob<KnobsTF1::RangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
}

FitLandauPol0::~FitLandauPol0()
{
}

void FitLandauPol0::Draw()
{
    func->Draw("same");
}

void FitLandauPol0::Fit(TH1 *hist)
{
    FitFunction::doFit(hist);
}

void FitLandauPol0::SetDefaults(TH1 *hist)
{
    func->SetParameter(3, 0); // offset
    func->SetParameter(2, 3.0);

    if(hist) {
        func->SetParameter(0, hist->GetMaximum()/3.0);
        const double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
        func->SetParameter(1,max_pos);

    } else {
        func->SetParameter(0,1000);
        func->SetParameter(1,0);
    }
}

void FitLandauPol0::SetRange(ant::interval<double> i)
{
    setRange(func, i);
}

ant::interval<double> FitLandauPol0::GetRange() const
{
    return getRange(func);
}

FitFunction::SavedState_t FitLandauPol0::Save() const
{
    std::vector<double> params;

    saveTF1(func,params);

    return params;
}

void FitLandauPol0::Load(const SavedState_t &data)
{
    if(data.size() != std::size_t(2+func->GetNpar())) {
        LOG(WARNING) << "Can't load parameters";
        return;
    }
    SavedState_t::const_iterator pos = data.begin();
    loadTF1(pos, func);

    Sync();

}

double FitLandauPol0::GetPeakPosition() const
{
    return MPV_trafo(func->GetParameter(1), func);
}

double FitLandauPol0::GetPeakWidth() const
{
    return func->GetParameter(2);
}
