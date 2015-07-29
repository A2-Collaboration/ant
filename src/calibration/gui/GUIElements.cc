#include "GUIElements.h"

#include "TLatex.h"


using namespace ant::calibration::gui;


class UndeleteableTLatex: public TLatex {
public:
    using TLatex::TLatex;

    /* remove the delete option from the context menu */
    virtual void Delete(Option_t*) override {}
};




void IndicatorLine::update_me()
{
    SetPosition(knob.get());
}

IndicatorLine::IndicatorLine(ant::calibration::gui::VirtualKnob &k): TLine(),
    knob(k),
    label(new UndeleteableTLatex(fX1,fY1,knob.name.c_str()))
{
    SetBit(kCanDelete, false);
    SetLineWidth(3);
    SetColor(k.color);
}

IndicatorLine::~IndicatorLine()
{
}

void IndicatorLine::Delete(Option_t *)
{
}

void IndicatorLine::Draw(Option_t *option)
{
    updateLabel();
    TLine::Draw(option);
    label->Draw();
}

void IndicatorLine::SetColor(const Color_t color)
{
    SetLineColor(color);
    label->SetTextColor(color);
}


void VerticalIndicatorLine::updateLabel()
{
    label->SetX(GetX1());
    label->SetY(GetY1());
    label->SetTextAngle(90);
}

VerticalIndicatorLine::VerticalIndicatorLine(VirtualKnob &k):
    IndicatorLine(k)
{
    SetPosition(k.get());
    updateLabel();
}

void VerticalIndicatorLine::SetX1(Double_t x1)
{
    fX1=x1;fX2=x1; updateLabel(); knob.set(x1); update_other();
}

void VerticalIndicatorLine::SetX2(Double_t x2) {fX2=x2;fX1=x2; updateLabel(); knob.set(x2); update_other();}

void VerticalIndicatorLine::SetY1(Double_t) {updateLabel();}

void VerticalIndicatorLine::SetY2(Double_t) {updateLabel();}

void VerticalIndicatorLine::SetupY(Double_t y1, Double_t y2) {fY1 = y1; fY2=y2; updateLabel(); }

void VerticalIndicatorLine::SetPosition(double p) { fX1=p; fX2=p; updateLabel(); }

double VerticalIndicatorLine::GetPosition() const { return GetX1(); }

void VerticalIndicatorLine::RangeUpdate(const Viewport &p) {
    SetupY(p.y1,p.y2);
}


void HorizontalIndicatorLine::updateLabel() {
    label->SetX(GetX1());
    label->SetY(GetY1());
    label->SetTextAngle(0);
}

HorizontalIndicatorLine::HorizontalIndicatorLine(VirtualKnob &k): IndicatorLine(k) { SetPosition(k.get()); updateLabel(); }

void HorizontalIndicatorLine::SetY1(Double_t y1) {fY1=y1;fY2=y1; updateLabel(); knob.set(y1); update_other();}

void HorizontalIndicatorLine::SetY2(Double_t y2) {fY2=y2;fY1=y2; updateLabel(); knob.set(y2); update_other();}

void HorizontalIndicatorLine::SetX1(Double_t) {updateLabel();}

void HorizontalIndicatorLine::SetX2(Double_t) {updateLabel();}

void HorizontalIndicatorLine::SetupX(Double_t x1, Double_t x2) {fX1 = x1; fX2=x2; updateLabel(); }

void HorizontalIndicatorLine::SetPosition(double p) { fY1=p; fY2=p; updateLabel(); }

double HorizontalIndicatorLine::GetPosition() const { return GetY1(); }

void HorizontalIndicatorLine::RangeUpdate(const Viewport &p) {
    SetupX(p.x1,p.x2);
}
