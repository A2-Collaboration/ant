#pragma once

#include "GUIbase.h"

#include <string>

class TF1;

namespace ant {
namespace calibration {
namespace gui {
namespace KnobsTF1 {

class ParameterKnob: public VirtualKnob {
protected:
    TF1* func = nullptr;
    const int parameter_index = 0;
public:

    ParameterKnob(const std::string& n, TF1* Func, int par, GUIElementDescription::GUI_Type type, Color_t color=kBlue, double LineW=3);

    virtual double get() const override;
    virtual void set(double a) override;

};

class ReferenceParameterKnob: public VirtualKnob {
protected:
    TF1* func = nullptr;
    const int parameter_index = 0;
    const int ref_index = 0;
public:

    ReferenceParameterKnob(const std::string& Name, TF1* Func, int par, int reference, GUIElementDescription::GUI_Type type, Color_t color=kBlue, double LineW=3);

    virtual double get() const override;
    virtual void set(double a) override;
    virtual double reference() const override;

};

class RangeKnob: public VirtualKnob {
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

}
}
}
}
