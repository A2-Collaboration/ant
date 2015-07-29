#include "FitFunction.h"

#include "base/interval.h"

#include "TF1.h"
#include "TH1.h"

using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;




FitFunctionGaus::FitFunctionGaus(double A, double x0, double sigma, interval<double> range):
    f(new TF1("","gaus", range.Start(), range.Stop()))
{
    f->SetParameter(0,A);
    f->SetParameter(1,x0);
    f->SetParameter(2,sigma);
    f->SetNpx(1000);

    knobs.emplace_back(&knob_A);
    knobs.emplace_back(&knob_x0);
    knobs.emplace_back(&knob_w);
}

FitFunctionGaus::~FitFunctionGaus()
{
}

void FitFunctionGaus::Draw()
{
    f->Draw("same");
}

FitFunction::knoblist_t FitFunctionGaus::getKnobs()
{
    return knobs;
}

void FitFunctionGaus::Fit(TH1 *hist)
{
    hist->Fit(f,"RBQN");
}


FitFunctionGaus::MyKnob::MyKnob(const std::string &n, TF1 *func, int par, VirtualKnob::GUI_Type gui):
    VirtualKnob(n,gui),
    f(func),
    p(par)
{
}

double FitFunctionGaus::MyKnob::get() const
{
    return f->GetParameter(p);
}

void FitFunctionGaus::MyKnob::set(double a)
{
    f->SetParameter(p,a);
}


FitFunctionGaus::MyWKnob::MyWKnob(const std::string &n, TF1 *func):
    VirtualKnob(n,GUI_Type::slider_vertical),
    f(func)
{
}

double FitFunctionGaus::MyWKnob::get() const
{
    return f->GetParameter(1) + f->GetParameter(2);
}

void FitFunctionGaus::MyWKnob::set(double a)
{
    auto v = a - f->GetParameter(1);
    f->SetParameter(2,v);
}
