#pragma once

#include "calibration/gui/GUIbase.h"

#include <memory>

#include "TLine.h"

class TLatex;

namespace ant {
namespace calibration {
namespace gui {

class IndicatorLine: public TLine, public GUIIndicator {
protected:
    ant::calibration::gui::VirtualKnob& knob;

    std::unique_ptr<TLatex> label;

    virtual void updateLabel() =0;

    virtual void update_other() {
        if(update)
            update->UpdateMe();
    }

public:
    IndicatorLine(VirtualKnob& k);
    virtual ~IndicatorLine();

    virtual void UpdateMe() override;

    update_notify_traits* update = nullptr;



    virtual void Delete(Option_t*) override;
    virtual void Draw(Option_t *option="") override;

    virtual void SetColor(const Color_t color);

};

class VerticalIndicatorLine: public IndicatorLine {
public:
protected:
    virtual void updateLabel() override;

public:
    VerticalIndicatorLine(VirtualKnob& k);

    virtual ~VerticalIndicatorLine() = default;

    // move both points the same way in x
    virtual void SetX1(Double_t x1) override;
    virtual void SetX2(Double_t x2) override;

    // ignore all y movements
    virtual void SetY1(Double_t) override;
    virtual void SetY2(Double_t) override;

    // set y positions
    virtual void SetupY(Double_t y1, Double_t y2);

    virtual void SetPosition(double p) override;
    virtual double GetPosition() const override;

    virtual void RangeUpdate(const Viewport& p) override;


};

class HorizontalIndicatorLine: public IndicatorLine {
public:
protected:
    virtual void updateLabel() override;

public:
    HorizontalIndicatorLine(VirtualKnob& k);

    virtual ~HorizontalIndicatorLine() = default;

    // move both points the same way in y
    virtual void SetY1(Double_t y1) override;
    virtual void SetY2(Double_t y2) override;

    // ignore all x movements
    virtual void SetX1(Double_t) override;
    virtual void SetX2(Double_t) override;

    // set x positions
    virtual void SetupX(Double_t x1, Double_t x2);

    virtual void SetPosition(double p) override;
    virtual double GetPosition() const override;

    virtual void RangeUpdate(const Viewport& p) override;

};


}
}
}
