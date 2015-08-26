#include "ManagerWindow.h"

#include "Manager.h"
#include "CalCanvas.h"

#include "TExec.h"
#include "TRootEmbeddedCanvas.h"
#include "TApplication.h"

#include "TGNumberEntry.h"
#include "TGStatusBar.h"
#include "TGButton.h"
#include "TGProgressBar.h"

#include "base/Logger.h"

#include <type_traits>

using namespace std;

namespace ant {
namespace calibration {
namespace gui {

class EmbeddedCanvas : public TRootEmbeddedCanvas {
public:
    EmbeddedCanvas(const TGWindow *p = 0) :
        TRootEmbeddedCanvas(0, p, 200, 200)
    {
        auto frame = (TGCompositeFrame*)fCanvasContainer;
        frame->RemoveInput(kKeyPressMask | kKeyReleaseMask);
    }
};

template<class theWidget>
class ActionWidget : public theWidget {
    struct MyExec : TExec {
        MyExec(function<void()> action_) : action(action_) {}
        virtual void Exec(const char*) override {
            action();
        }
    private:
        function<void()> action;
    };
    unique_ptr<MyExec> exec;
    bool* ptr_flag = nullptr;

public:
    using theWidget::theWidget;

    void SetAction(function<void()> action) {
        if(exec)
            return;
        exec = std_ext::make_unique<MyExec>(action);
        TQObject::Connect("Clicked()", "TExec", exec.get(), "Exec(=\"\")");
    }

    void LinkFlag(bool& flag) {
        static_assert(is_same<theWidget, TGCheckButton>::value, "LinkFlag only makes sense for check buttons");
        ptr_flag = addressof(flag);
        SetAction([this] () {
            *this->ptr_flag = this->IsOn();
        });
        SetFlag(flag);
    }

    void SetFlag(bool flag) {
        if(ptr_flag == nullptr)
            return;
        this->SetState(flag ? EButtonState::kButtonDown : EButtonState::kButtonUp);
        *this->ptr_flag = flag;
    }
};

void ManagerWindow::CreateToolbar(TGVerticalFrame* frame)
{
    // first row  with loop control commands

    TGHorizontalFrame* frm1 = new TGHorizontalFrame(frame,200,40);

    auto btn_autocontinue = new ActionWidget<TGCheckButton>(frm1,"AutoContinue");
    btn_autocontinue->LinkFlag(Mode.autoContinue);

    auto btn_showfit = new ActionWidget<TGCheckButton>(frm1,"Show each fit");
    btn_showfit->LinkFlag(Mode.showEachFit);

    auto btn_prev = new ActionWidget<TGTextButton>(frm1,"Prev (b)");
    keys[kKey_b] = btn_prev;
    btn_prev->SetAction([this, btn_autocontinue] () {
        btn_autocontinue->SetFlag(false);
        Mode.channelStep = -1;
        Mode.gotoNextSlice = false;
        manager->Run();
    });

    auto btn_next = new ActionWidget<TGTextButton>(frm1,"Next (n)");
    keys[kKey_n] = btn_next;
    btn_next->SetAction([this, btn_autocontinue] () {
        btn_autocontinue->SetFlag(false);
        Mode.channelStep = 1;
        Mode.gotoNextSlice = false;
        manager->Run();
    });

    auto btn_goto = new ActionWidget<TGTextButton>(frm1,"Goto");

    TGNumberEntry* numberentry = new TGNumberEntry(frm1, 0, 3, -1,
                                                   TGNumberFormat::kNESInteger,
                                                   TGNumberFormat::kNEANonNegative
                                                   );

    auto btn_finish = new ActionWidget<TGTextButton>(frm1,"Finish Slice");
    btn_finish->SetAction([this, btn_autocontinue] () {
        btn_autocontinue->SetFlag(true);
        Mode.channelStep = 1;
        Mode.gotoNextSlice = true;
        manager->Run();
    });

    // second row with fit specific commands
    /// \todo Make those commands specific for one canvas if
    /// more than one fitfunction is displayed...?!

    TGHorizontalFrame* frm2 = new TGHorizontalFrame(frame,200,40);

    auto btn_fit = new ActionWidget<TGTextButton>(frm2,"Fit (f)");
    keys[kKey_f] = btn_fit;
    btn_fit->SetAction([this] () {
        for(auto canvas : canvases)
            canvas->Fit();
    });

    auto btn_defaults = new ActionWidget<TGTextButton>(frm2,"SetDefaults (d)");
    keys[kKey_d] = btn_defaults;
    btn_defaults->SetAction([this] () {
        for(auto canvas : canvases)
            canvas->SetDefaults();
    });

    auto btn_undopop = new ActionWidget<TGTextButton>(frm2,"Undo pop (u)");
    keys[kKey_u] = btn_undopop;
    btn_undopop->SetAction([this] () {
        for(auto canvas : canvases)
            canvas->UndoPop();
    });

    auto btn_undopush = new ActionWidget<TGTextButton>(frm2,"Undo push (i)");
    keys[kKey_i] = btn_undopush;
    btn_undopush->SetAction([this] () {
        for(auto canvas : canvases)
            canvas->UndoPush();
    });

    // add them all together...

    auto layout_btn = new TGLayoutHints(kLHintsLeft,2,2,2,2);
    frm1->AddFrame(btn_prev, layout_btn);
    frm1->AddFrame(btn_next, layout_btn);
    frm1->AddFrame(btn_goto, layout_btn);
    frm1->AddFrame(numberentry, layout_btn);
    frm1->AddFrame(btn_finish, layout_btn);
    frm1->AddFrame(btn_autocontinue, layout_btn);
    frm1->AddFrame(btn_showfit, layout_btn);


    frm2->AddFrame(btn_fit, layout_btn);
    frm2->AddFrame(btn_defaults, layout_btn);
    frm2->AddFrame(btn_undopop, layout_btn);
    frm2->AddFrame(btn_undopush, layout_btn);


    // some progress bars

    TGHProgressBar* progress1 = new TGHProgressBar(frame);
    TGHProgressBar* progress2 = new TGHProgressBar(frame);

    auto layout_frm =  new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0);
    frame->AddFrame(frm1, layout_frm);
    frame->AddFrame(frm2, layout_frm);
    frame->AddFrame(progress1, layout_frm);
    frame->AddFrame(progress2, layout_frm);
}

void ManagerWindow::UpdateLayout()
{
    // Map all subwindows of main frame
    MapSubwindows();
    Resize(GetDefaultSize()); // this is used here to init layout algorithm
    MapWindow();
}

ManagerWindow::ManagerWindow(Manager* manager_) :
    TGMainFrame(gClient->GetRoot(), 400, 400),
    manager(manager_)
{
    // Set a name to the main frame
    SetWindowName("Ant-calib GUI");

    TGVerticalFrame* frame = new TGVerticalFrame(this);

    // Create a horizontal frame widget with buttons
    CreateToolbar(frame);

    // Create frame for canvases
    frame_canvases = new TGHorizontalFrame(frame,200,40);
    frame->AddFrame(frame_canvases, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX, 0, 0, 0, 0));

    // Statusbar
    Int_t parts[] = {45, 15, 10, 30};
    statusbar = new TGStatusBar(frame, 50, 10, kVerticalFrame);
    statusbar->SetParts(parts, 4);
    statusbar->Draw3DCorner(kFALSE);
    frame->AddFrame(statusbar, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 10, 0));

    AddFrame(frame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

    // after everthing is setup,
    manager->InitCanvases(this);

    AddInput(kKeyPressMask | kKeyReleaseMask);
    UpdateLayout();

    // set focus
    gVirtualX->SetInputFocus(GetId());
}

Bool_t ManagerWindow::HandleKey(Event_t* event) {

    if (event->fType == kGKeyPress) {
        char input[10];
        UInt_t keysym;
        gVirtualX->LookupString(event, input, sizeof(input), keysym);

        auto it_key = keys.find((EKeySym)keysym);
        if(it_key != keys.end()) {
            it_key->second->Clicked();
            return kTRUE;
        }
    }
    return TGMainFrame::HandleKey(event);
}

CalCanvas* ManagerWindow::AddCalCanvas(const string& name) {
    auto ecanvas = new EmbeddedCanvas(frame_canvases);
    auto canvas = new gui::CalCanvas(name.c_str(),ecanvas->GetCanvasWindowId());
    canvas->ConnectStatusBar(statusbar);
    ecanvas->AdoptCanvas(canvas);
    frame_canvases->AddFrame(ecanvas, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX, 0,0,0,0));
    canvases.push_back(canvas);
    UpdateLayout();
    return canvas;
}

ManagerWindow::~ManagerWindow()
{
    gApplication->Terminate(0);
}



}}} // namespace ant::calibration::gui

