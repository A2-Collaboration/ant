#include "EditorCanvas.h"

#include "calibration/Editor.h"
#include "calibration/gui/EditorWindow.h"
#include "base/std_ext/string.h"

#include "TH1D.h"
#include "TROOT.h"


#include <iostream>

using namespace ant;
using namespace std;
using namespace ant::calibration::gui;

EmbeddedEditorCanvas::EmbeddedEditorCanvas(EditorWindow* EditorWindow, const TGWindow* p) :
    TRootEmbeddedCanvas(0, p, 400, 400) // only important place to set some width/height
{
    auto frame = (TGCompositeFrame*)fCanvasContainer;
    frame->RemoveInput(kKeyPressMask | kKeyReleaseMask);
    theCanvas = new EditorCanvas(EditorWindow, GetCanvasWindowId());
    AdoptCanvas(theCanvas);
}

void EmbeddedEditorCanvas::EditSelection()
{
    theCanvas->StartEditData();
}

void EmbeddedEditorCanvas::SetToAverage()
{
    theCanvas->SetToAverage();
}

void EmbeddedEditorCanvas::ResetData()
{
    theCanvas->ResetCalibration();
}

void EmbeddedEditorCanvas::ApplyChanges()
{
    theCanvas->ApplyDataChanges();
}

void EmbeddedEditorCanvas::UpdateMe()
{
    theCanvas->UpdateMe();
}




EditorCanvas::EditorCanvas(EditorWindow* EditorWindow, int winID ):
    TCanvas("Editor",10,10,winID),
    editorWindow(EditorWindow),
    editor(EditorWindow->GetEditor())
{}

void EditorCanvas::UpdateMe()
{
    Modified();
    Update();
}



void EditorCanvas::ResetCalibration()
{
    editor->ResetData();
    StartEditData();
}

void EditorCanvas::ApplyDataChanges()
{
    for (auto i = 0; i < calDataHist->GetNbinsX() ; ++i)
        editor->cdata.Data.at(i).Value  = calDataHist->GetBinContent(i+1);
    editor->Save();
}

void EditorCanvas::StartEditData()
{
    gROOT->SetEditHistograms(kTRUE);
    calDataHist = new TH1D( (std_ext::formatter() << "hist-" << editor->cdata.CalibrationID
                             ).str().c_str(),
                            (std_ext::formatter() << "Data for " << editor->cdata.CalibrationID
                             ).str().c_str(),
                        editor->cdata.Data.size(), 0, editor->cdata.Data.size()
                        );

    for (const auto& entry: editor->cdata.Data)
        calDataHist->SetBinContent(entry.Key+1,entry.Value);

    calDataHist->SetStats(false);
    calDataHist->SetMarkerStyle(20);
    calDataHist->Draw("P");

    UpdateMe();
}

void EditorCanvas::HandleInput(EEventType button, Int_t x, Int_t y)
{
    TCanvas::HandleInput(button,x,y);
}

void EditorCanvas::SetToAverage()
{
    double mean = 0;
    for ( int i = 1; i <= calDataHist->GetNbinsX(); ++i)
        mean += calDataHist->GetBinContent(i);

    mean = 1.0 * mean / calDataHist->GetNbinsX();

    for ( int i = 1; i <= calDataHist->GetNbinsX(); ++i)
        calDataHist->SetBinContent(i,mean);

    UpdateMe();
}
