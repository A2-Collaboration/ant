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

struct GUIElementDescription {

    enum class GUI_Type {
        slider_vertical,
        slider_horizontal,
        textbox
    };

    GUIElementDescription(GUI_Type type, Color_t color, double lineWidth): Type(type), Color(color), LineWidth(lineWidth) {}

    GUI_Type Type;
    Color_t Color;
    double LineWidth;
};

class VirtualKnob {
public:

    std::string name;

    GUIElementDescription gui;


    VirtualKnob(const std::string& name_, GUIElementDescription Gui):
        name(name_),
        gui(Gui)
    {}

    // no copies of knobs!
    VirtualKnob(const VirtualKnob&) = delete;
    VirtualKnob& operator=(const VirtualKnob&) = delete;
    VirtualKnob(VirtualKnob&&) = delete;
    VirtualKnob& operator=(VirtualKnob&&) = delete;

    virtual double get() const =0;
    virtual void set(double v) =0;
};

}
}
}
