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
    virtual void update_me() =0;
};


class GUIIndicator: public update_notify_traits {
public:
    virtual ~GUIIndicator() = default;

    virtual void SetPosition(double p) =0;
    virtual double GetPosition() const =0;

    virtual void RangeUpdate(const Viewport& p) =0;
};


class VirtualKnob {
public:
    enum class GUI_Type {
        slider_vertical,
        slider_horizontal,
        textbox
    };

    std::string name;
    GUI_Type GUI;
    Color_t color;

    VirtualKnob(const std::string& name_="", const GUI_Type gui=GUI_Type::textbox, const Color_t col=kBlue):
        name(name_),
        GUI(gui),
        color(col)
    {}

    virtual double get() const =0;
    virtual void set(double v) =0;
};

}
}
}
