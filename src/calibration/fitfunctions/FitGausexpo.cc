#include "FitGausexpo.h"

#include "BaseFunctions.h"


#include "base/TF1Ext.h"

#include "TF1.h"
#include "TH1.h"

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


    func->SetParName(0,"A");
    func->SetParLimits(0, 0.0, 1E+12);
    func->SetParName(1, "x_{0}");
    func->SetParName(2, "sigma");
    func->SetParLimits(2, 0.0, 1E+12);
    func->SetParName(3, "p_{0}");
    func->SetParName(4, "p_{1}");
    func->SetParName(5, "p_{2}");




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
    FitFunction::doFit(hist, func);
    Sync();
}

void gui::FitGausexpo::FitBackground(TH1* hist)
{
    const auto fixedPars ={0,1,2};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist,func);
    Sync();
    UnFixParameters(func, fixedPars);
}

void gui::FitGausexpo::FitSignal(TH1* hist)
{
    const auto fixedPars ={3,4,5};
    FixParameters(func, fixedPars);
    FitFunction::doFit(hist,func);
    Sync();
    UnFixParameters(func, fixedPars);

}

void ant::calibration::gui::FitGausexpo::SetDefaults(TH1 *hist)
{
    func->SetParameter(0, hist->GetMaximum());
    const double max_pos= hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());

    auto range = GetRange();
    func->SetParameter(1, range.Clip(max_pos));

    func->SetParameter(2, 12);
    func->SetParameter(3, 1);
    func->SetParameter(4, 1);
    func->SetParameter(5,1);

    Sync();
}

void ant::calibration::gui::FitGausexpo::SetRange(ant::interval<double> i)
{
    setRange(func,i);
    setRange(signal, i);
    setRange(bg, i);

    func ->SetParLimits(1, i.Start(),i.Stop());
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


void ant::calibration::gui::FitGausexpo::Load(const std::vector<double> &data)
{

}

double ant::calibration::gui::FitGausexpo::GetPeakPosition() const
{
    return func->GetParameter(1);
}

double ant::calibration::gui::FitGausexpo::GetPeakWidth() const
{
    return func->GetParameter(2);
}









