#include "EditorWindow.h"

#include "base/std_ext/memory.h"

#include "TExec.h"
#include "TApplication.h"
#include "TSystem.h"

#include "TGNumberEntry.h"
#include "TGStatusBar.h"
#include "TGButton.h"
#include "TROOT.h"

#include "base/Logger.h"

#include <type_traits>
#include <sstream>
#include <iomanip>

using namespace std;

namespace ant {
namespace calibration {
namespace gui {

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

void EditorWindow::CreateToolbar(TGVerticalFrame* frame)
{
    // first row  with loop control commands

    TGHorizontalFrame* frm1 = new TGHorizontalFrame(frame);

    auto btn_autocontinue = new ActionWidget<TGCheckButton>(frm1,"AutoContinue");
    btn_autocontinue->LinkFlag(flag);

    auto btn_prev = new ActionWidget<TGTextButton>(frm1,"Prev (b)");
    keys[kKey_b] = btn_prev;
    btn_prev->SetAction([this] () {
        cout << "prev" << endl;
    });

    calibSelector = new MyComboBox(frm1);
    calibSelector->SetList(editor.GetListOfCalibrations());
    auto btn_select = new ActionWidget<TGTextButton>(frm1,"Select");
    btn_select->SetAction([this] () {
        cout <<  this->calibSelector->GetSelectedText() << endl;
    });

    auto entry_gotochannel = new TGNumberEntry(frm1, 0, 3, -1,
                                               TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative
                                               );

    auto btn_goto = new ActionWidget<TGTextButton>(frm1,"Goto");
    btn_goto->SetAction([this, entry_gotochannel, btn_autocontinue] () {
        cout <<  entry_gotochannel->GetIntNumber() << endl;
    });

    // second row with fit specific commands
    /// \todo Make those commands specific for one canvas if
    /// more than one fitfunction is displayed...?!


    // add them all together...
    auto layout_btn = new TGLayoutHints(kLHintsLeft,2,2,2,2);

    auto add_to_frame = [this, layout_btn] (TGHorizontalFrame* frm, TGWidget* widget) {
        frm->AddFrame(dynamic_cast<TGFrame*>(widget), layout_btn);
    };

    add_to_frame(frm1, btn_prev);
    add_to_frame(frm1, btn_select);
    add_to_frame(frm1, calibSelector);
    add_to_frame(frm1, btn_goto);
    add_to_frame(frm1, entry_gotochannel);
    frm1->AddFrame(btn_autocontinue, layout_btn);


    auto layout_frm =  new TGLayoutHints(kLHintsTop | kLHintsExpandX);
    frame->AddFrame(frm1, layout_frm);
}

void EditorWindow::UpdateLayout()
{
    // Map all subwindows of main frame
    MapSubwindows();
    Resize(GetDefaultSize()); // this is used here to init layout algorithm
    MapWindow();
}

EditorWindow::EditorWindow(const string& folder) :
    TGMainFrame(gClient->GetRoot())
{
    // Set a name to the main frame
    SetWindowName( (std_ext::formatter() << "Ant-calib Editor: " << folder).str().c_str() );

    editor.AddFromFolder(folder);

    TGVerticalFrame* frame = new TGVerticalFrame(this);

    // Create a horizontal frame widget with buttons
    CreateToolbar(frame);

    AddFrame(frame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));


    AddInput(kKeyPressMask | kKeyReleaseMask);
    UpdateLayout();

    // set focus
    gVirtualX->SetInputFocus(GetId());

}

Bool_t EditorWindow::HandleKey(Event_t* event) {

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

EditorWindow::~EditorWindow()
{
    // executed if window is closed
    gApplication->Terminate(0);
}

void EditorWindow::EditorCanvas::UpdateMe()
{
    Modified();
    Update();
}

EditorWindow::MyComboBox::MyComboBox(const TGWindow* p, Int_t id, UInt_t options, Pixel_t back):
    TGComboBox(p,id,options,back){}

void EditorWindow::MyComboBox::SetList(const list<string>& items)
{
    int i = 0;
    for (const auto& item: items)
        AddEntry(item.c_str(),i++);
}

string EditorWindow::MyComboBox::GetSelectedText()
{
    TGTextLBEntry* tgl = (TGTextLBEntry*) GetSelectedEntry();
    return string(tgl->GetText()->GetString());
}



}}} // namespace ant::calibration::gui

