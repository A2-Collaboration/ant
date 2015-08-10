#include "CalCanvas.h"

#include "FitFunction.h"
#include "GUIElements.h"
#include "base/Logger.h"

#include "TF1.h"
#include "TH1.h"
#include "TRootCanvas.h"

using namespace ant::calibration::gui;


CalCanvas::CalCanvas(const std::string &name):
    TCanvas(name.c_str()) {

    rootcanvas = dynamic_cast<TRootCanvas*>(GetCanvasImp());
    rootcanvas->DontCallClose();
}

CalCanvas::~CalCanvas() {
    ClearInidators();
}

void CalCanvas::Show(TH1 *h, FitFunction* f) {

    // empty UndoStack
    while(!UndoStack.empty()) {
        UndoStack.pop();
    }

    func = f;
    func->SetPoints(1000);
    hist = h;

    this->cd();
    h->Draw();
    f->Draw();

    SetupGUI();

    Modified();
    Update();
}

void CalCanvas::UpdateMe() {
    for(auto& i : indicators) {
        i->UpdateMe();
    }
    func->Sync();
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
    switch (knob.gui.Type) {
    case GUIElementDescription::GUI_Type::slider_horizontal:
        return MakeHorizontalIndicatorLine(knob);
        break;
    case GUIElementDescription::GUI_Type::slider_vertical:
        return MakeVerticalIndicatorLine(knob);
        break;
    case GUIElementDescription::GUI_Type::textbox:
        //@todo Implement
        return MakeVerticalIndicatorLine(knob);
    }
    return MakeVerticalIndicatorLine(knob);
}

void CalCanvas::HandleKeypress(const char key)
{
    switch (key) {
    case 'f':
        Fit();
        break;
    case 'u':
        UndoPop();
        break;
    case 'i':
        UndoPush();
        break;
    case 'd':
        SetDefaults();
        break;
    default:
        break;
    }

    if(gui_mode==nullptr)
        return;

    switch (key) {
    case '\r':
        gui_mode->channelStep = 1;
        gui_mode->alwaysDisplayFit = false;
        gui_mode->gotoNextRange = true;
        rootcanvas->Emit("CloseWindow()");
        break;
    case 'n':
        gui_mode->channelStep = 1;
        gui_mode->alwaysDisplayFit = true;
        gui_mode->gotoNextRange = false;
        rootcanvas->Emit("CloseWindow()");
        break;
    case 'b': // brevious item haha...
        gui_mode->channelStep = -1;
        gui_mode->alwaysDisplayFit = true;
        gui_mode->gotoNextRange = false;
        rootcanvas->Emit("CloseWindow()");
        break;
    default:
        break;
    }

}

void CalCanvas::SetDefaults()
{
    if(func && hist) {
        UndoPush();
        func->SetDefaults(hist);
        Modified();
        Update();
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

        UndoPush();

        func->Fit(hist);

        Modified();
        Update();
    }
}

void CalCanvas::UndoPush()
{
    if(func) {
        VLOG(7) << "Saving state to undo stack";
        UndoStack.push(func->Save());
    }
}

void CalCanvas::UndoPop()
{
    if(func) {
        if(!UndoStack.empty()) {

            VLOG(7) << "Loading state from undo stack";

            func->Load(UndoStack.top());

            if(UndoStack.size()>1) {
                UndoStack.pop();
            }

            Modified();
            UpdateMe();
        } else {
            VLOG(7) << "No earlier states on the stack";
        }
    }
}

void CalCanvas::ConnectReturnFunc(const char* receiver_class, void* receiver, const char* slot)
{
    rootcanvas->Connect("CloseWindow()", receiver_class, receiver, slot);
}

void CalCanvas::LinkGUIMode(CalCanvasMode* guimode_)
{
    gui_mode = guimode_;
}

void CalCanvas::Execute(const char *method, const char *params, Int_t *error) {
    std::string cmd(method);
    if(cmd == "Fit") {
        Fit();
    }
    else
    {
        TCanvas::Execute(method,params,error);
    }
}

void CalCanvas::Update() {

    auto p = getViewport();
    for(auto& i : indicators) {
        i->RangeUpdate(p);
        i->UpdateMe();
    }

    TCanvas::Update();

}

void CalCanvas::HandleInput(EEventType button, Int_t x, Int_t y)
{
    TCanvas::HandleInput(button,x,y);

    if(button == kKeyPress)
        HandleKeypress(x);
}


