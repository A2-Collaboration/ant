
#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGaus.h"

#include "base/cbtaps_display/TH2CB.h"
#include "base/std_ext.h"


#include "TCanvas.h"
#include "TRint.h"
#include "TGClient.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TRandom.h"
#include "TGButton.h"
#include "TGFrame.h"
#include "TRootEmbeddedCanvas.h"
#include "KeySymbols.h"
#include "TExec.h"
#include "TH1D.h"
#include "TGStatusBar.h"

#include <iostream>
#include <functional>
#include <memory>

using namespace std;
using namespace ant;
using namespace ant::calibration;

class MyTextButton : public TGTextButton {

    struct MyExec : TExec {
        MyExec(std::function<void()> action_) : action(action_) {}
        virtual void Exec(const char*) override {
            action();
        }
    private:
        std::function<void()> action;
    };
    std::unique_ptr<MyExec> exec;

public:
    MyTextButton(const TGWindow *p, const std::string& label) :
        TGTextButton(p, label.c_str())
    {
    }
    void SetAction(std::function<void()> action) {
        if(exec)
            return;
        exec = std_ext::make_unique<MyExec>(action);
        Connect("Clicked()", "TExec", exec.get(), "Exec(=\"\")");
    }
};

class MyMainFrame : public TGMainFrame
{
private:
    std::list<gui::CalCanvas*> canvases;
    TGHorizontalFrame* frame_canvases;
    TGStatusBar* statusbar;
public:
    MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
    virtual Bool_t HandleKey(Event_t *event) override;
    virtual ~MyMainFrame();
    gui::CalCanvas* AddCalCanvas(const std::string& name = "");
};

Bool_t MyMainFrame::HandleKey(Event_t *event) {

    if (event->fType == kGKeyPress) {
        char input[10];
        UInt_t keysym;

        gVirtualX->LookupString(event, input, sizeof(input), keysym);

        switch ((EKeySym)keysym) {
        case kKey_Left:
            cout << "pressed left arrow" << endl;
            return kTRUE;
        default:
            break;
        }
    }
    return TGMainFrame::HandleKey(event);
}

MyMainFrame::~MyMainFrame()
{
    gApplication->Terminate(0);
}

gui::CalCanvas* MyMainFrame::AddCalCanvas(const string& name) {
    auto ecanvas = new TRootEmbeddedCanvas(0,frame_canvases,200,200);
    auto canvas = new gui::CalCanvas(name.c_str(),ecanvas->GetCanvasWindowId());
    canvas->ConnectStatusBar(statusbar);
    ecanvas->AdoptCanvas(canvas);
    frame_canvases->AddFrame(ecanvas, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX, 0,0,0,0));
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
    return canvas;
}

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h) :
    TGMainFrame(p, w, h)
{
    BindKey(GetParent(), gVirtualX->KeysymToKeycode(kKey_Left),0);
    AddInput(kKeyPressMask | kKeyReleaseMask);

    TGVerticalFrame* frame = new TGVerticalFrame(this);

    // Create a horizontal frame widget with buttons
    TGHorizontalFrame* frame_buttons = new TGHorizontalFrame(frame,200,40);
    MyTextButton* draw = new MyTextButton(frame_buttons,"&Draw");
    draw->SetAction([this] () {
        auto canvas = AddCalCanvas();
        canvas->cd();
        auto h = new TH1D("","h",100,0,300);
        auto f = new gui::FitGaus();
        f->SetDefaults(nullptr); // using empty h is not meaningful for testing
        canvas->Show(h, f);

        auto canvas_display = AddCalCanvas();
        canvas_display->cd();
        auto cb = new TH2CB("","CB");
        cb->Draw();
    });

    auto button_layout = new TGLayoutHints(kLHintsLeft,5,5,3,4);

    frame_buttons->AddFrame(draw, button_layout);
    TGTextButton* exit = new TGTextButton(frame_buttons,"&Exit",
                                          "gApplication->Terminate(0)");
    frame_buttons->AddFrame(exit, button_layout);

    frame->AddFrame(frame_buttons, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));

    // Create frame for canvases
    frame_canvases = new TGHorizontalFrame(frame,200,40);
    frame->AddFrame(frame_canvases, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX, 0, 0, 0, 0));

    // Statusbar
    // status bar
    Int_t parts[] = {45, 15, 10, 30};
    statusbar = new TGStatusBar(frame, 50, 10, kVerticalFrame);
    statusbar->SetParts(parts, 4);
    statusbar->Draw3DCorner(kFALSE);
    frame->AddFrame(statusbar, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 10, 0));


    AddFrame(frame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

    // Set a name to the main frame
    SetWindowName("Simple Example");

    // Map all subwindows of main frame
    MapSubwindows();

    Resize(GetDefaultSize()); // this is used here to init layout algorithm

    // Map main frame
    MapWindow();
    gVirtualX->SetInputFocus(GetId());
}

int main(int argc, char** argv)
{
    TRint app("guitest",&argc,argv);
    // MainFram is destroyed on close!
    new MyMainFrame(gClient->GetRoot(),400,400);
    app.Run();
}
