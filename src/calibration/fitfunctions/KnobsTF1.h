#pragma once

#include "calibration/gui/Indicator_traits.h"

#include <string>
#include <functional>

class TF1;

namespace ant {
namespace calibration {
namespace gui {
namespace KnobsTF1 {


class ParameterKnob : public IndicatorKnob {
protected:
    TF1* func = nullptr;
    const int parameter_index = 0;
public:

    ParameterKnob(const std::string& Name, TF1* Func, int par,
                  IndicatorProperties::Type_t type, Color_t color=kBlue, double LineW=3);

    virtual double get() const override;
    virtual void set(double a) override;

};

class FixedParameterKnob : public IndicatorKnob {
protected:
    TF1* func = nullptr;
    const int parameter_index = 0;
public:

    FixedParameterKnob(const std::string& Name, TF1* Func, int par,
                       IndicatorProperties::Type_t type, Color_t color=kAzure-4, double LineW=3);

    virtual double get() const override;
    virtual void set(double a) override;

};

class RangedParameterKnob : public ParameterKnob {
public:
    enum class ConstraintType {
        lowerThanMin,
        lowerThanMax,
        higherThanMin,
        higherThanMax
    };

protected:
    const ConstraintType constraint_type;
public:
    RangedParameterKnob(
            const std::string& Name, TF1* Func, int par,
            ConstraintType constraint_type_,
            IndicatorProperties::Type_t gui_type, Color_t color=kBlue, double LineW=3);

    virtual void set(double a) override;
};

class TransformedParameterKnob : public ParameterKnob {
public:
    using transformation_t = std::function<double(double, TF1*)>;

    TransformedParameterKnob(const std::string& Name, TF1* Func, int par,
                             transformation_t trafo,
                             transformation_t trafo_inverse,
                             IndicatorProperties::Type_t type, Color_t color=kBlue, double LineW=3);

    virtual double get() const override;
    virtual void set(double a) override;
protected:
    const transformation_t transformation;
    const transformation_t transformation_inverse;
};

class ReferenceParameterKnob : public TransformedParameterKnob {

public:

    ReferenceParameterKnob(const std::string& Name, TF1* Func, int par, int reference,
                           IndicatorProperties::Type_t type, Color_t color=kBlue, double LineW=3);
};

class RangeKnob: public IndicatorKnob {
public:
    enum class RangeEndType {
        upper,
        lower
    };

protected:
    TF1* func = nullptr;
    const RangeEndType type;
public:

    RangeKnob(const std::string& Name, TF1* Func, RangeEndType Type, Color_t color=kBlack, double LineW=1);

    virtual double get() const override;
    virtual void set(double a) override;

};

class FixedRangeKnob: public RangeKnob {
public:

    FixedRangeKnob(const std::string& Name, TF1* Func, RangeEndType Type, Color_t color=kGray+1, double LineW=1);

    virtual double get() const override;
    virtual void set(double a) override;

};

}
}
}
}
