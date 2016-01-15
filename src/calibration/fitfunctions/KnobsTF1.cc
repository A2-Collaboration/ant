#include "KnobsTF1.h"

#include "TF1.h"

using namespace ant::calibration::gui;
using namespace ant::calibration::gui::KnobsTF1;


ParameterKnob::ParameterKnob(const std::string& name, TF1* Func, int par, IndicatorProperties::Type_t type, Color_t color, double LineW):
    IndicatorKnob(name,IndicatorProperties(type,color,LineW)),
    func(Func),
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





RangeKnob::RangeKnob(const std::string& Name, TF1* Func, RangeEndType Type, Color_t color, double LineW)
    :IndicatorKnob(Name,{IndicatorProperties::Type_t::slider_vertical,color,LineW}),
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

ReferenceParameterKnob::ReferenceParameterKnob(const std::string& Name, TF1* Func,
                                               int par, int reference,
                                               IndicatorProperties::Type_t type, Color_t color, double LineW):
    ParameterKnob(Name,Func,par,type,color,LineW),
    ref_index(reference)
{
}

double ReferenceParameterKnob::get() const
{
    return func->GetParameter(parameter_index) + reference();
}

void ReferenceParameterKnob::set(double a)
{
    func->SetParameter(parameter_index, a-reference());
}

double ReferenceParameterKnob::reference() const
{
    return func->GetParameter(ref_index);
}


RangedParameterKnob::RangedParameterKnob(const std::string& n, TF1* Func, int par,
                                         RangedParameterKnob::ConstraintType constraint_type_,
                                         IndicatorProperties::Type_t gui_type, Color_t color, double LineW) :
ParameterKnob(n, Func, par, gui_type, color, LineW),
  constraint_type(constraint_type_)
{
}

void RangedParameterKnob::set(double a)
{
    double min, max;
    func->GetRange(min, max);

    switch(constraint_type) {
    case ConstraintType::lowerThanMin:
        if(min<a)
            return;
        break;
    case ConstraintType::lowerThanMax:
        if(max<a)
            return;
        break;
    case ConstraintType::higherThanMax:
        if(max>a)
            return;
        break;
    case ConstraintType::higherThanMin:
        if(min>a)
            return;
        break;
    default:
        return;
    }
    ParameterKnob::set(a);
}
