#include "FitFunction.h"

#include "base/interval.h"
#include "TF1Knobs.h"

#include "TF1.h"
#include "TH1.h"

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;




FitFunctionGaus::FitFunctionGaus(double A, double x0, double sigma, interval<double> range):
    func(new TF1("","gaus", range.Start(), range.Stop()))
{
    func->SetParameter(0,A);
    func->SetParameter(1,x0);
    func->SetParameter(2,sigma);
    func->SetNpx(1000);

    Addknob<KnobsTF1::ParameterKnob>("A", func, 0, VirtualKnob::GUI_Type::slider_horizontal);
    Addknob<KnobsTF1::ParameterKnob>("x_0",func,1);
    Addknob<MyWKnob>("#sigma",func);
    Addknob<KnobsTF1::RangeKnob>("Min", func, KnobsTF1::RangeKnob::RangeEndType::lower);
    Addknob<KnobsTF1::RangeKnob>("Max", func, KnobsTF1::RangeKnob::RangeEndType::upper);
}

FitFunctionGaus::~FitFunctionGaus()
{
}

void FitFunctionGaus::Draw()
{
    func->Draw("same");
}

void FitFunctionGaus::Fit(TH1 *hist)
{
    hist->Fit(func,"RBQN");
}


FitFunctionGaus::MyWKnob::MyWKnob(const std::string &n, TF1 *Func):
    VirtualKnob(n,GUI_Type::slider_vertical),
    func(Func)
{
}

double FitFunctionGaus::MyWKnob::get() const
{
    return func->GetParameter(1) + func->GetParameter(2);
}

void FitFunctionGaus::MyWKnob::set(double a)
{
    auto v = a - func->GetParameter(1);
    func->SetParameter(2,v);
}


FitFunction::~FitFunction()
{

}
