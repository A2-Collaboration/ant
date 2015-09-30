#include "EditorCanvas.h"

#include "calibration/gui/EditorWindow.h"
#include "calibration/Editor.h"
#include "base/std_ext/string.h"
#include "tree/TCalibrationData.h"

#include "TH2D.h"
#include "TH2.h"
#include "TROOT.h"

#include <list>
#include <cmath>

#include <iostream>
#include "TStyle.h"

using namespace ant;
using namespace std;
using namespace ant::calibration::gui;

EmbeddedEditorCanvas::EmbeddedEditorCanvas(EditorWindow* EditorWindow, const string& calID, const TGWindow* p) :
    TRootEmbeddedCanvas(0, p, 400, 400) // only important place to set some width/height
{
    auto frame = (TGCompositeFrame*)fCanvasContainer;
    frame->RemoveInput(kKeyPressMask | kKeyReleaseMask);
    theCanvas = new EditorCanvas(EditorWindow, calID, GetCanvasWindowId());
    AdoptCanvas(theCanvas);
}

void EmbeddedEditorCanvas::SelectInvalid()
{
    theCanvas->MarkInvalid();
}

void EmbeddedEditorCanvas::SetCalID(const string& calID)
{
    theCanvas->SetCalID(calID);
}

std::list<uint32_t> EmbeddedEditorCanvas::GetSelected() const
{
    return theCanvas->CreateSelectionList();
}

void EmbeddedEditorCanvas::clearSelections()
{
    theCanvas->ResetCalibration();
}

void EmbeddedEditorCanvas::EditSelection()
{
    theCanvas->EditData();
}

bool EmbeddedEditorCanvas::InDataEditMode() const
{
    return theCanvas->getDataEditorFlag();
}

void EmbeddedEditorCanvas::UpdateMe()
{
    theCanvas->UpdateMe();

}




EditorCanvas::EditorCanvas(EditorWindow* EditorWindow, const string& calID, int winID ):
    TCanvas("Editor",10,10,winID),
    editorWindow(EditorWindow),
    editor(EditorWindow->GetEditor()),
    flag_intervalStart_set(false),
    flag_data_editor(false)
{


    Int_t palette[5];
    palette[0] = 0;
    palette[1] = 32;
    palette[2] = 38;
    palette[3] = 46;
    palette[4] = 42;
    gStyle->SetPalette(5,palette);

    SetCalID(calID);

}

void EditorCanvas::SetCalID(const string& calID)
{
    currentCalID = calID;

    auto nSteps = editor->GetNumberOfSteps(currentCalID);

    calHist = new TH2D( (std_ext::formatter() << "hist-" << currentCalID).str().c_str(),
                        (std_ext::formatter() << "History for " << currentCalID).str().c_str(),
                        100, 0, 100,
                        nSteps,
                        0, nSteps
                        );
    calHist->SetXTitle("TID [%]");
    calHist->SetYTitle("Calibration Step");
    if ( nSteps < 50 )
        calHist->GetYaxis()->SetNdivisions(nSteps);
    calHist->Draw("col");
    calHist->SetStats(false);
    calHist->GetZaxis()->SetRangeUser(0,4);

    ResetCalibration();
}

bool EditorCanvas::getDataEditorFlag() const
{
    return flag_data_editor;
}

list<uint32_t> EditorCanvas::CreateSelectionList()
{
    list<uint32_t> theList;
    for ( const auto theIndex: indexMemory)
        theList.emplace_back(theIndex);
    return theList;
}

void EditorCanvas::UpdateMe()
{
    Modified();
    Update();
}

void EditorCanvas::updateCalHist()
{
    calHist->Reset();
    for (const auto& ran: editor->GetAllRanges(currentCalID))
        for (int i = floor(ran.second.Start()); i < ceil(ran.second.Stop()) ; ++i)
            calHist->SetBinContent(i+1,ran.first+1,2);

    for (const auto& ran: editor->GetAllValidRanges(currentCalID))
        for (int i = floor(ran.second.Start()); i < ceil(ran.second.Stop()) ; ++i)
            calHist->SetBinContent(i+1,ran.first+1,1);

    calHist->Draw("col");
    UpdateMe();
}


void EditorCanvas::ResetCalibration()
{
    flag_intervalStart_set = false;
    flag_data_editor = false;
    gROOT->SetEditHistograms(kFALSE);
    indexMemory.clear();
    updateCalHist();
}



void EditorCanvas::fillLine(uint32_t lineNumber)
{
    for (Int_t i = 0; i < 100; ++i )
        calHist->Fill(i,lineNumber,2);
}

void EditorCanvas::unFillLine(uint32_t lineNumber)
{
    for (Int_t i = 0; i < 100; ++i )
        calHist->Fill(i,lineNumber,-2);
}



void EditorCanvas::MarkInvalid()
{
    auto valids = editor->GetAllValidRanges(currentCalID);
    set<uint32_t> validset;
    for (const auto v: valids)
        validset.emplace(v.first);

    for (uint32_t i = 0 ; i < editor->GetNumberOfSteps(currentCalID) ; ++i)
        if ( validset.find(i) == validset.end() && indexMemory.find(i) == indexMemory.end() )
        {
            cout << "added invalid line " << i << endl;
            indexMemory.emplace(i);
            fillLine(i);
        }
    UpdateMe();
}
void EditorCanvas::applyDataChanges(TCalibrationData& theData)
{
    for (auto i = 0; i < calDataHist->GetNbinsX() ; ++i)
        theData.Data.at(i).Value  = calDataHist->GetBinContent(i+1);
    ResetCalibration();
}

void EditorCanvas::StartEditData(TCalibrationData& theData, uint32_t stepIndex)
{
    gROOT->SetEditHistograms(kTRUE);
    calDataHist = new TH1D( (std_ext::formatter() << "hist-" << currentCalID
                                                  << "-"<< stepIndex
                             ).str().c_str(),
                            (std_ext::formatter() << "Data for " << currentCalID
                                                  << ", step " <<  stepIndex
                             ).str().c_str(),
                        theData.Data.size(), 0, theData.Data.size()
                        );

    for (const auto& entry: theData.Data)
        calDataHist->SetBinContent(entry.Key+1,entry.Value);

    calDataHist->SetStats(false);
    calDataHist->SetMarkerStyle(20);
    calDataHist->Draw("P");

    UpdateMe();
}

bool EditorCanvas::EditData()
{
    if (indexMemory.size() != 1)
        return false;
    uint32_t stepIndex = *(indexMemory.begin());
    TCalibrationData&  theData = editor->ModifyItem(currentCalID,stepIndex);

    if (flag_data_editor)
    {
        applyDataChanges(theData);
        return true;
    }

    flag_data_editor = true;
    StartEditData(theData, stepIndex);
    return true;
}

void EditorCanvas::markInterval(Int_t y)
{
    TH2* h = dynamic_cast<TH2*>(fClickSelected);
    if (h){
        auto index = floor(AbsPixeltoY(y));

        if (!flag_intervalStart_set)
        {
            indexInterVal.Start() = index;
            fillLine(index);
            flag_intervalStart_set = true;
        }
        else
        {
            indexInterVal.Stop() = index;
            indexInterVal.MakeSane();
            updateCalHist();
            for (auto inInt = indexInterVal.Start(); inInt <= indexInterVal.Stop(); ++inInt)
            {
                indexMemory.emplace(inInt);
                fillLine(inInt);
            }
        }
    }
    UpdateMe();
}

void EditorCanvas::markLine(Int_t y)
{
    TH2* h = dynamic_cast<TH2*>(fClickSelected);
    if (h){
        auto step = floor(AbsPixeltoY(y));
        auto foundat = indexMemory.find(step);
        if ( foundat == indexMemory.end() )
        {
            indexMemory.emplace(step);
            fillLine(step);
        }
        else
        {
            indexMemory.erase(foundat);
            unFillLine(step);
        }
    }
    UpdateMe();
    editorWindow->UpdateMe();
}

void EditorCanvas::HandleInput(EEventType button, Int_t x, Int_t y)
{
    TCanvas::HandleInput(button,x,y);

    if ( !flag_data_editor )
    {
        if (button == kButton1Down)
            markLine(y);
        //still buggy
//        if (button == kButton3Down)
//            markInterval(y);
    }
}
