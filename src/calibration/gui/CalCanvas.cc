#include "CalCanvas.h"

#include "Indicator.h"
#include "base/Logger.h"

#include "fitfunctions/FitFunction.h"

#include "TF1.h"
#include "TH1.h"
#include "TRootCanvas.h"
#include "TGStatusBar.h"

#include <sstream>

using namespace std;
using namespace ant::calibration::gui;


CalCanvas::CalCanvas(const std::string &name, const std::string& title):
    TCanvas(name.c_str(), title.c_str())
{

    rootcanvas = dynamic_cast<TRootCanvas*>(GetCanvasImp());
    rootcanvas->DontCallClose();
}

CalCanvas::CalCanvas(const string& name, Int_t winid) :
    TCanvas(name.c_str(), 10, 10, winid)
{

}

CalCanvas::~CalCanvas() {
    ClearIndicators();
}

void CalCanvas::Show(TH1* h, FitFunction* f, bool preserveYaxis) {

    // empty UndoStack
    while(!UndoStack.empty()) {
        UndoStack.pop();
    }

    this->cd();
    h->Draw();
    f->Draw();

    // preserve zoom state of axis
    if(hist != nullptr) {
        PreserveAxis(hist->GetXaxis(), h->GetXaxis());
        if(preserveYaxis) {
            // we cannot use PreserveAxis, since the Y axis is unbinned
            // for a one-dimensional histogram...
            Double_t ymin = hist->GetMinimum();
            Double_t ymax = hist->GetMaximum();
            h->SetAxisRange(ymin, ymax, "Y");
        }
    }

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

void CalCanvas::PreserveAxis(TAxis* axis1, TAxis* axis2)
{
    Int_t binmin = axis1->GetFirst();
    Int_t binmax = axis1->GetLast();
    Float_t xmin = axis1->GetBinLowEdge(binmin);
    Float_t xmax = axis1->GetBinLowEdge(binmax);
    Int_t newmin = axis2->FindBin(xmin);
    Int_t newmax = axis2->FindBin(xmax);
    axis2->SetRange(newmin,newmax);
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
}

void CalCanvas::Fit() {
    if(func && hist) {
        VLOG(3) << "Refitting";

        UndoPush();

        func->Fit(hist);

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
    TCanvas::Update();

    auto p = getViewport();
    for(auto& i : indicators) {
        i->RangeUpdate(p);
        i->UpdateMe();
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


