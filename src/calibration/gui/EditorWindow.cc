#include "EditorWindow.h"

#include "base/std_ext/memory.h"
#include "calibration/gui/EditorCanvas.h"

#include "TExec.h"
#include "TApplication.h"
#include "TSystem.h"

#include "TGNumberEntry.h"
#include "TGStatusBar.h"
#include "TGButton.h"
#include "TROOT.h"

#include "base/Logger.h"
#include "TH2D.h"

#include <type_traits>
#include <sstream>
#include <iomanip>
#include <TCanvas.h>

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

void EditorWindow::createSelector(TGVerticalFrame* frame)
{
    TGHorizontalFrame* frm1 = new TGHorizontalFrame(frame);

    calibSelector = new MyComboBox(frm1,0);
    calibSelector->SetList(editor->GetListOfCalibrations());
    auto btn_select = new ActionWidget<TGTextButton>(frm1,"Select");
    btn_select->SetAction([this] () {
        currentCalID = this->calibSelector->GetSelectedText();
        ecanvas->SetCalID(currentCalID);
    });

    // add them all together...
    auto layout_btn = new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsExpandY,2,2,2,2);

    auto add_to_frame = [this, layout_btn] (TGHorizontalFrame* frm, TGWidget* widget) {
        frm->AddFrame(dynamic_cast<TGFrame*>(widget), layout_btn);
    };

    add_to_frame(frm1, calibSelector);
    add_to_frame(frm1, btn_select);

    auto layout_frm =  new TGLayoutHints(kLHintsTop | kLHintsExpandX);
    frame->AddFrame(frm1, layout_frm);
}


void EditorWindow::createToolbar(TGVerticalFrame* frame)
{
    TGHorizontalFrame* frm2 = new TGHorizontalFrame(frame);
    TGHorizontalFrame* frm3 = new TGHorizontalFrame(frame);
    TGHorizontalFrame* frm4 = new TGHorizontalFrame(frame);


    auto btn_selectInValid = new ActionWidget<TGTextButton>(frm2,"Select invalid");
    keys[kKey_s] = btn_selectInValid;
    rootButton_markInValid = btn_selectInValid;
    btn_selectInValid->SetAction([this] () {
        this->ecanvas->SelectInvalid();
        UpdateMe();
    });

    auto btn_expandSelection = new ActionWidget<TGTextButton>(frm2,"expand selection");
    rootButton_expandSelection = btn_expandSelection;
    btn_expandSelection->SetAction([this] ()
    {
        this->ecanvas->ExpandSelection();
        this->ecanvas->clearSelections();
        UpdateMe();
    });

    auto btn_clear = new ActionWidget<TGTextButton>(frm3,"Clear selection / go back");
    keys[kKey_c] = btn_clear;
    btn_clear->SetAction([this] () {
        this->ecanvas->clearSelections();
        UpdateMe();
    });

    auto btn_edit = new ActionWidget<TGTextButton>(frm3,"start Editor / write data");
    rootButton_StartEditor = btn_edit;
    btn_edit->SetAction([this] () {
        this->ecanvas->EditSelection();
    });

    auto btn_dublicateLast = new ActionWidget<TGTextButton>(frm2,"dublicate last step");
    rootButton_dublicateLast = btn_dublicateLast;
    btn_dublicateLast->SetAction([this] () {
        this->editor->DublicateLast(currentCalID);
        this->ecanvas->SetCalID(currentCalID);
        UpdateMe();
    });

    auto btn_avg = new ActionWidget<TGTextButton>(frm3,"AVG");
    rootButton_avg = btn_avg;
    btn_avg->SetAction([this] () {
        this->ecanvas->SetToAverage();
    });

    auto btn_delete = new ActionWidget<TGTextButton>(frm2,"delete selection");
    rootButton_delete = btn_delete;
    btn_delete->SetAction([this] () {
        this->deleteSelections();
        UpdateMe();
    });

    auto btn_save = new ActionWidget<TGTextButton>(frm4,"save database");
    btn_save->SetAction([this] () {
        this->editor->SaveToFolder(dataFolder);
    });

    auto btn_Quit = new ActionWidget<TGTextButton>(frm4,"Exit without saving");
    btn_Quit->SetAction([this] () {
        gApplication->Terminate(0);
    });

    auto btn_saveQuit = new ActionWidget<TGTextButton>(frm4,"Save and Exit");
    btn_saveQuit->SetAction([this] () {
        this->editor->SaveToFolder(dataFolder);
        gApplication->Terminate(0);
    });


    // add them all together...
    auto layout_btn = new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsExpandY,2,2,2,2);
    auto add_to_frame = [this, layout_btn] (TGHorizontalFrame* frm, TGWidget* widget) {
        frm->AddFrame(dynamic_cast<TGFrame*>(widget), layout_btn);
    };

    add_to_frame(frm2, btn_selectInValid);
    add_to_frame(frm3, btn_edit);
    add_to_frame(frm2, btn_expandSelection);
    add_to_frame(frm2, btn_delete);
    add_to_frame(frm2, btn_dublicateLast);
    add_to_frame(frm3, btn_avg);
    add_to_frame(frm3, btn_clear);
    add_to_frame(frm4, btn_save);
    add_to_frame(frm4, btn_Quit);
    add_to_frame(frm4, btn_saveQuit);

    auto layout_frm =  new TGLayoutHints(kLHintsBottom | kLHintsExpandX);
    frame->AddFrame(frm4, layout_frm);
    frame->AddFrame(frm3, layout_frm);
    frame->AddFrame(frm2, layout_frm);
}

void EditorWindow::updateLayout()
{
    // Map all subwindows of main frame
    MapSubwindows();
    Resize(GetDefaultSize()); // this is used here to init layout algorithm
    MapWindow();
}

void EditorWindow::deleteSelections()
{
    for (const auto& i: ecanvas->GetSelected())
        editor->Remove(currentCalID,i);
    // has to load current Calibration as new one, because size changed:
    ecanvas->SetCalID(currentCalID);
}

void EditorWindow::disableButtons()
{
    rootButton_markInValid->SetEnabled(kTRUE);
    if (ecanvas->InDataEditMode())
        rootButton_markInValid->SetEnabled(kFALSE);

    rootButton_expandSelection->SetEnabled(ecanvas->GetSelected().size() > 0);

    rootButton_delete->SetEnabled(kTRUE);
    if ( ecanvas->GetSelected().size() == 0 || ecanvas->InDataEditMode())
        rootButton_delete->SetEnabled(kFALSE);

    rootButton_StartEditor->SetEnabled(kFALSE);
    if ( ecanvas->GetSelected().size() == 1)
        rootButton_StartEditor->SetEnabled(kTRUE);

    rootButton_dublicateLast->SetEnabled( ecanvas->GetSelected().size() == 0);

    rootButton_avg->SetEnabled( ecanvas->InDataEditMode());

}

EditorWindow::EditorWindow(const string& folder) :
    TGMainFrame(gClient->GetRoot()),
    editor(make_shared<ant::calibration::Editor>()),
    dataFolder(folder)
{
    // Set a name to the main frame
    SetWindowName( (std_ext::formatter() << "Ant-calib Editor: " << folder).str().c_str() );

    editor->AddFromFolder(folder);
    currentCalID = editor->GetListOfCalibrations().front();

    TGVerticalFrame* frame = new TGVerticalFrame(this);

    // Create a horizontal frame widget with buttons
    createSelector(frame);

    frame_canvas = new TGHorizontalFrame(frame);
    frame->AddFrame(frame_canvas, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX));

    ecanvas = new EmbeddedEditorCanvas(this,currentCalID,frame_canvas);
    frame_canvas->AddFrame(ecanvas, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX));
    updateLayout();

    AddFrame(frame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    createToolbar(frame);


    AddInput(kKeyPressMask | kKeyReleaseMask);
    UpdateMe();

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

std::shared_ptr<Editor> EditorWindow::GetEditor() { return editor;}


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

void EditorWindow::UpdateMe()
{
    disableButtons();
    updateLayout();
}



}}} // namespace ant::calibration::gui

