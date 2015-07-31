#include "TF1Knobs.h"

#include "TF1.h"

using namespace ant::calibration::gui;
using namespace ant::calibration::gui::KnobsTF1;


ParameterKnob::ParameterKnob(const std::string &n, TF1 *func, int par, VirtualKnob::GUI_Type gui):
    VirtualKnob(n,gui),
    func(func),
    parameter_index(par)
{
}

double ParameterKnob::get() const
{
    return func->GetParameter(parameter_index);
}

void ParameterKnob::set(double a)
{
    func->SetParameter(parameter_index,a);
}





RangeKnob::RangeKnob(const std::string& Name, TF1* Func, RangeEndType Type)
    :VirtualKnob(Name,GUI_Type::slider_vertical,kBlack),
      func(Func),
      type(Type)
{
}

double RangeKnob::RangeKnob::get() const
{
    double min, max;
    func->GetRange(min,max);

    if(type==RangeEndType::lower) {
        return min;
    }

    return max;
}

void RangeKnob::RangeKnob::set(double a)
{
    double min, max;
    func->GetRange(min,max);

    switch (type) {

    case RangeEndType::lower:
        if(a < max)
            func->SetRange(a,max);
        break;

    case RangeEndType::upper:
            if(a>min)
          func->SetRange(min,a);
        break;

    default:
        break;
    }
}
