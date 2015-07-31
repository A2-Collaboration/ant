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
    knobs.emplace_back(&knob_minR);
    knobs.emplace_back(&knob_maxR);
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


FitFunctionGaus::MyRangeKnob::MyRangeKnob(const std::string &n, TF1 *func, FitFunctionGaus::MyRangeKnob::RangeEndType t_)
    :VirtualKnob(n,GUI_Type::slider_vertical,kBlack),
      f(func),
      t(t_)
{
}

double FitFunctionGaus::MyRangeKnob::get() const
{
    double min, max;
    f->GetRange(min,max);

    if(t==RangeEndType::lower) {
        return min;
    }

    return max;
}

void FitFunctionGaus::MyRangeKnob::set(double a)
{
    double min, max;
    f->GetRange(min,max);

    switch (t) {

    case RangeEndType::lower:
        if(a < max)
            f->SetRange(a,max);
        break;

    case RangeEndType::upper:
            if(a>min)
          f->SetRange(min,a);
        break;

    default:
        break;
    }
}
