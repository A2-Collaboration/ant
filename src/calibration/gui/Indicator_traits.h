#pragma once

#include "Rtypes.h"

namespace ant {
namespace calibration {
namespace gui {

struct Viewport {
    Viewport(double x1_=.0, double y1_=.0, double x2_=.0, double y2_=.0):
        x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
    double x1;
    double y1;
    double x2;
    double y2;
};


class update_notify_traits {
public:
    virtual void UpdateMe() =0;
};


class Indicator: public update_notify_traits {
public:
    virtual ~Indicator() = default;

    virtual void SetPosition(double p) =0;
    virtual double GetPosition() const =0;

    virtual void RangeUpdate(const Viewport& p) =0;
};

struct IndicatorProperties {

    enum class Type_t {
        slider_vertical,
        slider_horizontal
    };

    IndicatorProperties(Type_t type, Color_t color, double lineWidth): Type(type), Color(color), LineWidth(lineWidth) {}

    Type_t Type;
    Color_t Color;
    double LineWidth;
};

class IndicatorKnob {
public:

    std::string name;
    IndicatorProperties gui;


    IndicatorKnob(const std::string& name_, IndicatorProperties Gui):
        name(name_),
        gui(Gui)
    {}

    // no copies of knobs!
    IndicatorKnob(const IndicatorKnob&) = delete;
    IndicatorKnob& operator=(const IndicatorKnob&) = delete;
    IndicatorKnob(IndicatorKnob&&) = delete;
    IndicatorKnob& operator=(IndicatorKnob&&) = delete;

    virtual double get() const =0;
    virtual void set(double v) =0;
    virtual double reference() const { return 0.0; }
};

}
}
}
