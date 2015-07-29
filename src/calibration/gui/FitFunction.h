#pragma once

#include "base/interval.h"
#include "calibration/gui/GUIbase.h"

#include <list>
#include <string>

#include "Rtypes.h"

class TH1;
class TF1;

namespace ant {
namespace calibration {
namespace gui {


class FitFunction {
public:
    using knoblist_t = std::list<ant::calibration::gui::VirtualKnob*>;

    virtual ~FitFunction() = default;
    virtual void Draw() =0;
    virtual knoblist_t getKnobs() =0;
    virtual void Fit(TH1* hist) =0;
};



class FitFunctionGaus: public FitFunction {
protected:
    TF1* f = nullptr;

    class MyKnob: public VirtualKnob {
    protected:
        TF1* f = nullptr;
        int p = 0;
    public:

        MyKnob(const std::string& n, TF1* func, int par, GUI_Type gui=GUI_Type::slider_vertical);

        virtual double get() const override;
        virtual void set(double a) override;

    };

    class MyWKnob: public VirtualKnob {
    protected:
        TF1* f = nullptr;
    public:

        MyWKnob(const std::string& n, TF1* func);

        virtual double get() const override;
        virtual void set(double a) override;

    };

    MyKnob knob_A  = MyKnob("A",f,0, VirtualKnob::GUI_Type::slider_horizontal);
    MyKnob knob_x0 = MyKnob("x_{0}",f,1);
    MyWKnob knob_w = MyWKnob("#sigma",f);

public:

    std::list<VirtualKnob*> knobs;

    FitFunctionGaus(double A, double x0, double sigma, ant::interval<double> range);

    virtual ~FitFunctionGaus();

    virtual void Draw() override;
    virtual knoblist_t getKnobs();

    virtual void Fit(TH1* hist) override;


};

}
}
}

