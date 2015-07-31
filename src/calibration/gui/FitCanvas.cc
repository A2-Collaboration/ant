#include "FitCanvas.h"

#include "FitFunction.h"
#include "GUIElements.h"
#include "base/Logger.h"

#include "TF1.h"
#include "TH1.h"


using namespace ant::calibration::gui;




CalCanvas::CalCanvas(const std::string &name):
    TCanvas(name.c_str()) {
}

CalCanvas::~CalCanvas() {
    ClearInidators();
}

void CalCanvas::Show(TH1 *h, std::shared_ptr<FitFunction> f) {
    func = f;
    hist = h;
    this->cd();
    h->Draw();
    f->Draw();
    Update();
    SetupGUI();
}

void CalCanvas::update_me() {
    for(auto& i : indicators) {
        i->update_me();
    }
    Update();
}

Viewport CalCanvas::getViewport() {
    Viewport p;
    GetRangeAxis(p.x1,p.y1,p.x2,p.y2);
    return p;
}

GUIIndicator *CalCanvas::MakeVerticalIndicatorLine(VirtualKnob &knob) {

    auto p = getViewport();

    auto tmp = new VerticalIndicatorLine(knob);
    tmp->update = this;
    tmp->RangeUpdate(p);
    tmp->Draw();

    return tmp;
}

GUIIndicator *CalCanvas::MakeHorizontalIndicatorLine(VirtualKnob &knob) {

    auto p = getViewport();

    auto tmp = new HorizontalIndicatorLine(knob);
    tmp->update = this;
    tmp->RangeUpdate(p);
    tmp->Draw();

    return tmp;
}

GUIIndicator *CalCanvas::MakeGUIElement(VirtualKnob &knob)
{
    switch (knob.GUI) {
    case VirtualKnob::GUI_Type::slider_horizontal:
        return MakeHorizontalIndicatorLine(knob);
        break;
    case VirtualKnob::GUI_Type::slider_vertical:
        return MakeVerticalIndicatorLine(knob);
        break;
    default:
        return MakeVerticalIndicatorLine(knob);
        break;
    }
}

void CalCanvas::HandleKeypress(const char key)
{
    switch (key) {
    case 'f':
        Fit();
        break;
    default:
        break;
    }
}

void CalCanvas::ClearInidators() {
    for(auto& i : indicators) {
        delete i;
    }
    indicators.clear();
}

void CalCanvas::SetupGUI() {

    ClearInidators();

    for(auto& knob : func->getKnobs()) {
        auto gui = MakeGUIElement(*knob);
        indicators.emplace_back(gui);
    }
}

/* remove snaping guides when moving objects on the canvas */
void CalCanvas::ShowGuidelines(TObject *, const Int_t, const char, const bool) {}


void CalCanvas::Fit() {
    if(func && hist) {
        VLOG(3) << "Refitting";

        func->Fit(hist);

        Modified();
        Update();
    }
}

void CalCanvas::Execute(const char *method, const char *params, Int_t *error) {
    std::string cmd(method);
    if(cmd == "Fit") {
        Fit();
    } else
    {
        TCanvas::Execute(method,params,error);
    }
}

void CalCanvas::Update() {

    auto p = getViewport();
    for(auto& i : indicators) {
        i->RangeUpdate(p);
        i->update_me();
    }

    TCanvas::Update();

}

void CalCanvas::HandleInput(EEventType button, Int_t x, Int_t y)
{

    if(button == kKeyPress) {
        HandleKeypress(x);
    } else {
        TCanvas::HandleInput(button,x,y);
    }
}
