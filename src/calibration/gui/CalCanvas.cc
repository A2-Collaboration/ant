#include "CalCanvas.h"

#include "Indicator.h"
#include "base/Logger.h"

#include "fitfunctions/FitFunction.h"

#include "TF1.h"
#include "TH1.h"
#include "TRootCanvas.h"
#include "TGStatusBar.h"
#include "TExec.h"

#include <sstream>

using namespace std;
using namespace ant::calibration::gui;



CalCanvas::~CalCanvas() {
    ClearIndicators();
}

void CalCanvas::Show(TH1* h, FitFunction* f, bool preserveYaxis) {

    // empty UndoStack
    while(!UndoStack.empty()) {
        UndoStack.pop();
    }

    // restore zoom state of axis
    if(axis_settings.HaveX()) {
        TAxis* xaxis = h->GetXaxis();
        xaxis->SetRange(xaxis->FindBin(axis_settings.x_min),
                        xaxis->FindBin(axis_settings.x_max));
    }
    if(preserveYaxis && axis_settings.HaveY()) {
        // we cannot use PreserveAxis, since the Y axis is unbinned
        // for a one-dimensional histogram...
        h->SetAxisRange(axis_settings.y_min, axis_settings.y_max, "Y");
    }

    this->cd();
    h->Draw();
    f->Draw();

    func = f;
    hist = h;

    SetupGUI();
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

Indicator* CalCanvas::MakeVerticalIndicatorLine(IndicatorKnob &knob) {

    auto p = getViewport();

    auto tmp = new VerticalIndicatorLine(knob);
    tmp->update = this;
    tmp->RangeUpdate(p);
    tmp->Draw();

    return tmp;
}

Indicator* CalCanvas::MakeHorizontalIndicatorLine(IndicatorKnob &knob) {

    auto p = getViewport();

    auto tmp = new HorizontalIndicatorLine(knob);
    tmp->update = this;
    tmp->RangeUpdate(p);
    tmp->Draw();

    return tmp;
}

Indicator* CalCanvas::MakeGUIElement(IndicatorKnob &knob)
{
    switch (knob.gui.Type) {
    case IndicatorProperties::Type_t::slider_horizontal:
        return MakeHorizontalIndicatorLine(knob);
        break;
    case IndicatorProperties::Type_t::slider_vertical:
        return MakeVerticalIndicatorLine(knob);
        break;
    }
    return nullptr;
}

void CalCanvas::ClearIndicators() {
    for(auto& i : indicators) {
        delete i;
    }
    indicators.clear();
}

void CalCanvas::SetupGUI() {

    // the indicators need an updated canvas
    // in order to figure out the right viewport
    Modified();
    Update();

    ClearIndicators();

    for(auto& knob : func->GetKnobs()) {
        auto gui = MakeGUIElement(*knob);
        indicators.emplace_back(gui);
    }

    // afterwards we update to actually show
    // the created indicators
    Modified();
    Update();

    class ExecUpdate : public TExec {
    protected:
        CalCanvas* canvas = nullptr;
    public:
        ExecUpdate(CalCanvas* c) : canvas(c)        {}
        virtual void Exec(const char*) override {
            canvas->Update(true);
        }
    };

    hist->GetListOfFunctions()->Add(new ExecUpdate(this));
}

void CalCanvas::Fit(const FitType_t type) {
    if(func && hist) {
        VLOG(3) << "Refitting";

        UndoPush();

        switch (type) {
            case FitType_t::Total:
                func->Fit(hist);
                break;
            case FitType_t::Signal:
                func->FitSignal(hist);
                break;
            case FitType_t::Background:
                func->FitBackground(hist);
                break;
        }

        Modified();
        Update();

        UpdateMe();

        Modified();
        Update();
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

void CalCanvas::Update() {
    Update(false);
}

void CalCanvas::Update(bool fromhist)
{
    TCanvas::Update();
    auto p = getViewport();
    for(auto& i : indicators) {
        i->RangeUpdate(p);
        i->UpdateMe();
    }

    // remember the axis settings
    if(fromhist && hist != nullptr) {
        TAxis* xaxis = hist->GetXaxis();
        axis_settings.x_min = xaxis->GetBinLowEdge(xaxis->GetFirst());
        axis_settings.x_max = xaxis->GetBinLowEdge(xaxis->GetLast());

        axis_settings.y_min = hist->GetMinimum();
        axis_settings.y_max = hist->GetMaximum();
    }
}

void CalCanvas::ConnectStatusBar(TGStatusBar* statusbar_)
{
    statusbar = statusbar_;
}

void CalCanvas::ProcessedEvent(Int_t event, Int_t x, Int_t y, TObject* selected)
{
    if(statusbar) {
        statusbar->SetText(selected->GetTitle(), 0);
        statusbar->SetText(selected->GetName(),  1);
        stringstream ss;
        if(event == kKeyPress)
           ss << (char)x;
        else
           ss << x << "," << y;
        statusbar->SetText(ss.str().c_str(), 2);
        statusbar->SetText(selected->GetObjectInfo(x,y), 3);
    }
    TCanvas::ProcessedEvent(event, x, y, selected);
}


