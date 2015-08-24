#include "base/std_ext.h"
#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGaus.h"

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

#include <iostream>
#include <functional>
#include <memory>

using namespace std;
using namespace ant;

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
    calibration::gui::CalCanvas* canvas;
public:
    MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
    virtual Bool_t HandleKey(Event_t *event) override;
    virtual ~MyMainFrame();
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

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h) :
    TGMainFrame(p, w, h)
{
    BindKey(GetParent(), gVirtualX->KeysymToKeycode(kKey_Left),0);
    AddInput(kKeyPressMask | kKeyReleaseMask);

    // Create a horizontal frame widget with buttons
    TGHorizontalFrame *hframe = new TGHorizontalFrame(this,200,40);
    MyTextButton *draw = new MyTextButton(hframe,"&Draw");
    draw->SetAction([this] () {
        canvas->cd();
        auto h = new TH1D("h","h",100,0,300);
        auto f = new calibration::gui::FitGaus();
        f->SetDefaults(nullptr); // using empty h is not meaningful for testing
        canvas->Show(h, f);
    });

    hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX,
                                             5,5,3,4));
    TGTextButton *exit = new TGTextButton(hframe,"&Exit",
                                          "gApplication->Terminate(0)");
    hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));

    AddFrame(hframe, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));

    // Create canvas widget
    auto ecanvas = new TRootEmbeddedCanvas(0,this,200,200);
    canvas = new calibration::gui::CalCanvas("calcanvas",ecanvas->GetCanvasWindowId());
    ecanvas->AdoptCanvas(canvas);
    AddFrame(ecanvas, new TGLayoutHints(kLHintsExpandX |
                                        kLHintsExpandY, 10,10,10,10));

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
