#include "ManagerWindow.h"

#include "Manager.h"
#include "CalCanvas.h"

#include "TExec.h"
#include "TRootEmbeddedCanvas.h"
#include "TApplication.h"
#include "TSystem.h"

#include "TGNumberEntry.h"
#include "TGStatusBar.h"
#include "TGButton.h"
#include "TGProgressBar.h"
#include "TROOT.h"

#include "base/Logger.h"

#include <type_traits>
#include <sstream>
#include <iomanip>

using namespace std;

namespace ant {
namespace calibration {
namespace gui {

class EmbeddedCanvas : public TRootEmbeddedCanvas {
public:
    EmbeddedCanvas(const TGWindow *p = 0) :
        TRootEmbeddedCanvas(0, p, 400, 400) // only important place to set some width/height
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

class ProgressBar : public TGHProgressBar {
    string label;
public:
    ProgressBar(const TGWindow* p, const string& label_) :
        TGHProgressBar(p),
        label(label_)
    {
        ShowPosition(true, false, label.c_str());
    }

    void SetValue(unsigned value) {
        Reset();  // reset fixes some superweird when going backwards...
        SetPosition(value);


        auto max = static_cast<unsigned>(GetMax());

        stringstream ss;
        if(value<=max) {
            stringstream ss_max;
            ss_max << max;
            ss << label << " "
               << setw(ss_max.str().size()) << setfill('0') << value
               << "/" << ss_max.str();
        }
        else {
            ss << "Finished " << label;
        }
        // we misuse the format string...
        ShowPosition(true, false, ss.str().c_str());
    }

};


void ManagerWindow::CreateToolbar(TGVerticalFrame* frame)
{
    // first row  with loop control commands

    TGHorizontalFrame* frm1 = new TGHorizontalFrame(frame);

    auto btn_autocontinue = new ActionWidget<TGCheckButton>(frm1,"AutoContinue");
    btn_autocontinue->LinkFlag(Mode.autoContinue);

    auto btn_showfit = new ActionWidget<TGCheckButton>(frm1,"Show each fit");
    btn_showfit->LinkFlag(Mode.showEachFit);

    auto btn_prev = new ActionWidget<TGTextButton>(frm1,"Prev (b)");
    keys[kKey_b] = btn_prev;
    btn_prev->SetAction([this] () {
        Mode.channelStep = -1;
        Mode.gotoNextSlice = false;
        RunManager();
    });

    auto btn_next = new ActionWidget<TGTextButton>(frm1,"Next (n)");
    keys[kKey_n] = btn_next;
    btn_next->SetAction([this] () {
        Mode.channelStep = 1;
        Mode.gotoNextSlice = false;
        RunManager();
    });


    auto entry_gotochannel = new TGNumberEntry(frm1, 0, 3, -1,
                                               TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative
                                               );

    auto btn_goto = new ActionWidget<TGTextButton>(frm1,"Goto");
    btn_goto->SetAction([this, entry_gotochannel, btn_autocontinue] () {
        btn_autocontinue->SetFlag(false);
        Mode.gotoNextSlice = false;
        Mode.requestChannel = entry_gotochannel->GetIntNumber();
        RunManager();
    });

    auto btn_finish = new ActionWidget<TGTextButton>(frm1,"Finish Slice (space)");
    keys[kKey_Space] = btn_finish;
    btn_finish->SetAction([this, btn_autocontinue] () {
        Mode.channelStep = 1;
        Mode.gotoNextSlice = true;
        RunManager();
    });

    // second row with fit specific commands
    /// \todo Make those commands specific for one canvas if
    /// more than one fitfunction is displayed...?!

    TGHorizontalFrame* frm2 = new TGHorizontalFrame(frame);

    auto btn_fit = new ActionWidget<TGTextButton>(frm2,"Fit (f)");
    keys[kKey_f] = btn_fit;
    btn_fit->SetAction([this] () {
        for(auto canvas : canvases)
            canvas->Fit();
    });

    auto btn_fitsignal = new ActionWidget<TGTextButton>(frm2,"Fit Signal (s)");
    keys[kKey_s] = btn_fitsignal;
    btn_fitsignal->SetAction([this] () {
        for(auto canvas : canvases)
            canvas->Fit(CalCanvas::FitType_t::Signal);
    });

    auto btn_fitbackground = new ActionWidget<TGTextButton>(frm2,"Fit Background (B)");
    keys[kKey_B] = btn_fitbackground;
    btn_fitbackground->SetAction([this] () {
        for(auto canvas : canvases)
            canvas->Fit(CalCanvas::FitType_t::Background);
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

    auto btn_nopreviousparams = new ActionWidget<TGCheckButton>(frm2,"Ignore previous params");
    btn_nopreviousparams->LinkFlag(Mode.FitOptions.IgnorePreviousFitParameters);

    // add them all together...
    auto layout_btn = new TGLayoutHints(kLHintsLeft,2,2,2,2);

    auto add_nonfinish = [this, layout_btn] (TGHorizontalFrame* frm, TGWidget* widget) {
        frm->AddFrame(dynamic_cast<TGFrame*>(widget), layout_btn);
        nonfinish_widgets.push_back(widget);
    };

    add_nonfinish(frm1, btn_prev);
    add_nonfinish(frm1, btn_next);
    add_nonfinish(frm1, btn_goto);
    add_nonfinish(frm1, entry_gotochannel);
    frm1->AddFrame(btn_finish, layout_btn);
    frm1->AddFrame(btn_autocontinue, layout_btn);
    frm1->AddFrame(btn_showfit, layout_btn);


    add_nonfinish(frm2, btn_fit);
    add_nonfinish(frm2, btn_fitsignal);
    add_nonfinish(frm2, btn_fitbackground);
    add_nonfinish(frm2, btn_defaults);
    add_nonfinish(frm2, btn_undopop);
    add_nonfinish(frm2, btn_undopush);
    add_nonfinish(frm2, btn_nopreviousparams);

    // some progress bars
    progress_channel = new ProgressBar(frame, "Channel");
    progress_slice = new ProgressBar(frame, "Slice");

    auto layout_frm =  new TGLayoutHints(kLHintsTop | kLHintsExpandX);
    frame->AddFrame(frm1, layout_frm);
    frame->AddFrame(frm2, layout_frm);
    frame->AddFrame(progress_channel, layout_frm);
    frame->AddFrame(progress_slice, layout_frm);


}

void ManagerWindow::UpdateLayout()
{
    // Map all subwindows of main frame
    MapSubwindows();
    Resize(GetDefaultSize()); // this is used here to init layout algorithm
    MapWindow();
}

void ManagerWindow::RunManager()
{
    if(running)
        return;
    running  = true;
    while(true) {
        auto ret = manager->Run();

        if(ret == Manager::RunReturn_t::Wait) {
            for(auto canvas : canvases)
                canvas->Update();
            if(!gROOT->IsBatch())
                break;
        }
        else if(ret == Manager::RunReturn_t::Exit) {
            if(!gROOT->IsBatch())
                gApplication->Terminate(0);
            break;
        }
        gSystem->ProcessEvents();
    }
    running = false;
}

ManagerWindow::ManagerWindow(Manager* manager_) :
    TGMainFrame(gClient->GetRoot()),
    manager(manager_)
{
    if(gROOT->IsBatch()) {
        Mode.channelStep = 1;
        Mode.gotoNextSlice = true;
        Mode.autoContinue = true;
        manager->InitGUI(this);
        RunManager();
        return;
    }

    // Set a name to the main frame
    SetWindowName("Ant-calib GUI");

    TGVerticalFrame* frame = new TGVerticalFrame(this);

    // Create a horizontal frame widget with buttons
    CreateToolbar(frame);

    // Create frame for canvases
    frame_canvases = new TGHorizontalFrame(frame);
    frame->AddFrame(frame_canvases, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX));

    // Statusbar
    Int_t parts[] = {45, 15, 10, 30};
    statusbar = new TGStatusBar(frame, 50, 10, kVerticalFrame);
    statusbar->SetParts(parts, 4);
    statusbar->Draw3DCorner(kFALSE);
    frame->AddFrame(statusbar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

    AddFrame(frame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));


    AddInput(kKeyPressMask | kKeyReleaseMask);
    UpdateLayout();

    // set focus
    gVirtualX->SetInputFocus(GetId());

    // after everthing is setup,
    // init the manager
    manager->InitGUI(this);
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
    gui::CalCanvas* canvas = nullptr;

    stringstream canvas_name;
    if(name.empty()) {
        canvas_name << "CalCanvas_" << canvases.size();
    }
    else {
        canvas_name << name;
    }

    if(gROOT->IsBatch()) {
        canvas = new gui::CalCanvas(canvas_name.str());
    }
    else {
        auto ecanvas = new EmbeddedCanvas(frame_canvases);
        canvas = new gui::CalCanvas(canvas_name.str(), ecanvas->GetCanvasWindowId());
        canvas->ConnectStatusBar(statusbar);
        ecanvas->AdoptCanvas(canvas);
        frame_canvases->AddFrame(ecanvas, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX));
        UpdateLayout();
    }
    canvases.push_back(canvas);
    return canvas;
}

void ManagerWindow::SetProgressMax(unsigned slices, unsigned channels)
{
    if(gROOT->IsBatch())
        return;
    progress_slice->SetRange(0, slices);
    progress_channel->SetRange(0, channels);
}

void ManagerWindow::SetProgress(unsigned slice, unsigned channel)
{
    if(gROOT->IsBatch())
        return;
    progress_slice->SetValue(slice);
    progress_channel->SetValue(channel);
}

void ManagerWindow::SetFinishMode(bool flag)
{
    for(auto w : nonfinish_widgets) {
        auto btn = dynamic_cast<TGButton*>(w);
        if(btn != nullptr) {
            btn->SetEnabled(!flag);
            continue;
        }
        auto entry = dynamic_cast<TGTextEntry*>(w);
        if(entry != nullptr) {
            entry->SetEnabled(!flag);
        }
    }
}

ManagerWindow::~ManagerWindow()
{
    // executed if window is closed
    gApplication->Terminate(0);
}



}}} // namespace ant::calibration::gui

