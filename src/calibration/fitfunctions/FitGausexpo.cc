#include "FitGausexpo.h"

#include "BaseFunctions.h"


#include "base/TF1Ext.h"

#include "TF1.h"
#include "TH1.h"
#include "TSpectrum.h"

#include <cmath>

using namespace std;
using namespace ant::calibration;

void ant::calibration::gui::FitGausexpo::Sync()
{
    signal->SetParameters(&(func->GetParameters()[0]));
    bg->SetParameters(    &(func->GetParameters()[3]));
    setRange(signal, GetRange());
    setRange(bg, GetRange());



}

ant::calibration::gui::FitGausexpo::FitGausexpo()
{
    signal= functions::gaus::getTF1();
    signal->SetLineColor(kRed);

    bg = functions::exponential::getTF1();
    bg->SetLineColor(kBlue);

    func = functions::Gausexpo::getTF1();
    func->SetLineColor(kGreen);




    SetRange(ant::interval<double>(25,140));
    func->SetParName(0,"A");
    func->SetParLimits(0, 0.0, 5000);
    func->SetParName(1, "x_{0}");
    func->SetParName(2, "sigma");
    func->SetParLimits(2, 0.0, 1E+12);
    func->SetParName(3, "p_{0}");
    func->SetParName(4, "p_{1}");




    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(0), func, 0, IndicatorProperties::Type_t::slider_horizontal);
    AddKnob<KnobsTF1::ParameterKnob>(func->GetParName(1), func, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::ReferenceParameterKnob>(func->GetParName(2), func, 2, 1, IndicatorProperties::Type_t::slider_vertical);
    AddKnob<KnobsTF1::RangeKnob>("min",func,KnobsTF1::RangeKnob::RangeEndType::lower);
    AddKnob<KnobsTF1::RangeKnob>("max",func,KnobsTF1::RangeKnob::RangeEndType::upper);


}

ant::calibration::gui::FitGausexpo::~FitGausexpo()
{
    delete signal;
    delete bg;
    delete func;


}

void ant::calibration::gui::FitGausexpo::Draw()
{
    signal->Draw("same");
    bg->Draw("same");
    func->Draw("same");
}
void ant::calibration::gui::FitGausexpo::Fit(TH1* hist)
{
    FitFunction::doFit(hist);
    Sync();
}

void gui::FitGausexpo::FitBackground(TH1* hist)
{
    const auto fixedPars ={0,1,2};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    Sync();
    UnFixParameters(func, fixedPars);
}

void gui::FitGausexpo::FitSignal(TH1* hist)
{
    const auto fixedPars ={3,4,5};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist);
    Sync();
    UnFixParameters(func, fixedPars);

}




void ant::calibration::gui::FitGausexpo::SetDefaults(TH1 *hist)
{
    /* Idee um Peak zu finden, noch nicht ideal...  */
//      double xp=65;
//      double yp=700;
//      TSpectrum* s = new TSpectrum();
//      double nfound = s->Search(hist);
//      float *xpeaks = s->GetPositionX();
//      for(int p=0; p<nfound; p++){
//             xp = xpeaks[p];
//             int bin = hist->GetXaxis()->FindBin(xp);
//             yp = hist->GetBinContent(bin);
//      }
//      if(xp<25 || xp>120)
//      {
//          xp=60;
//      }

    func->SetParameter(0, 700);
    const double max_pos= hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    auto range = GetRange();

    func->SetParameter(1,65);
    func ->SetParLimits(1, 20,100);

    func->SetParameter(2, 20);
    func->SetParameter(3, log( range.Clip(max_pos)));
    func->SetParameter(4, -0.05);


    Sync();
 //   delete s;
}

void ant::calibration::gui::FitGausexpo::SetRange(ant::interval<double> i)
{
    setRange(func,i);
    setRange(signal, i);
    setRange(bg, i);

//    func ->SetParLimits(1, 20,100);
//    func ->SetParLimits(0, 0,2000);
}

ant::interval<double> ant::calibration::gui::FitGausexpo::GetRange() const
{
    return getRange(func);
}

std::vector<double> ant::calibration::gui::FitGausexpo::Save() const
{
    SavedState_t params;
    params.reserve(2+func->GetNpar());
    saveTF1(func,params);

    return params;
}


void ant::calibration::gui::FitGausexpo::Load(const std::vector<double>&)
{

}

double ant::calibration::gui::FitGausexpo::GetPeakPosition() const
{
    return func->GetParameter(1);

}

double ant::calibration::gui::FitGausexpo::GetPeakError() const
{
    return func->GetParError(1);
}

double ant::calibration::gui::FitGausexpo::GetPeakWidth() const
{
    return func->GetParameter(2);
}

double ant::calibration::gui::FitGausexpo::GetPeakWidtherr() const
{
    return func->GetParError(2);
}









