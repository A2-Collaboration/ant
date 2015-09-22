#include "EditorCanvas.h"

#include "calibration/Editor.h"
#include "base/std_ext/string.h"

#include "TH2D.h"
#include "TH2.h"

#include <list>
#include <cmath>

#include <iostream>

using namespace ant;
using namespace std;
using namespace ant::calibration::gui;

EmbeddedEditorCanvas::EmbeddedEditorCanvas(const std::shared_ptr<calibration::Editor>& editor, const string& calID, const TGWindow* p) :
    TRootEmbeddedCanvas(0, p, 400, 400) // only important place to set some width/height
{
    auto frame = (TGCompositeFrame*)fCanvasContainer;
    frame->RemoveInput(kKeyPressMask | kKeyReleaseMask);
    theCanvas = new EditorCanvas(editor, calID, GetCanvasWindowId());
    AdoptCanvas(theCanvas);
}

void EmbeddedEditorCanvas::selectUnvalid()
{
    theCanvas->MarkUnValid();
}

void EmbeddedEditorCanvas::SetCalID(const string& calID)
{
    theCanvas->SetCalID(calID);
}

std::list<uint32_t> EmbeddedEditorCanvas::GetSelected()
{
    return theCanvas->CreateSelectionList();
}

void EmbeddedEditorCanvas::clearSelections()
{
    theCanvas->ResetCalibration();
}

void EmbeddedEditorCanvas::UpdateMe()
{
    theCanvas->UpdateMe();

}
EditorCanvas::EditorCanvas(const std::shared_ptr<Editor>& editor, const string& calID, int winID ):
    TCanvas("Editor",10,10,winID),
    ed(editor)
{
    SetCalID(calID);
}

void EditorCanvas::SetCalID(const string& calID)
{
    currentCalID = calID;

    //        if ( !calHist == nullptr )
    //           delete calHist;

    calHist = new TH2D( (std_ext::formatter() << "hist-" << currentCalID).str().c_str(),
                        (std_ext::formatter() << "History for " << currentCalID).str().c_str(),
                        100, 0, 100,
                        ed->GetNumberOfSteps(currentCalID),
                        0, ed->GetNumberOfSteps(currentCalID)
                        );

    calHist->SetXTitle("TID [%]");
    calHist->SetYTitle("Calibration Step");
    calHist->Draw("col");
    calHist->SetStats(false);
    updateCalHist();
    UpdateMe();
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
    for (const auto& ran: ed->GetAllRanges(currentCalID))
        for (int i = floor(ran.second.Start()); i < ceil(ran.second.Stop()) ; ++i)
            calHist->SetBinContent(i+1,ran.first+1,2);

    for (const auto& ran: ed->GetAllValidRanges(currentCalID))
        for (int i = floor(ran.second.Start()); i < ceil(ran.second.Stop()) ; ++i)
            calHist->SetBinContent(i+1,ran.first+1,1);

    UpdateMe();
}


void EditorCanvas::ResetCalibration()
{
    intervalStartSet = false;
    indexMemory.clear();
    updateCalHist();
}



void EditorCanvas::fillLine(uint32_t lineNumber)
{
    for (Int_t i = 0; i < 100; ++i )
        calHist->Fill(i,lineNumber,3);
}

void EditorCanvas::unFillLine(uint32_t lineNumber)
{
    for (Int_t i = 0; i < 100; ++i )
        calHist->Fill(i,lineNumber,-3);
}



void EditorCanvas::MarkUnValid()
{
    auto valids = ed->GetAllValidRanges(currentCalID);
    set<uint32_t> validset;
    for (const auto v: valids)
        validset.emplace(v.first);

    for (uint32_t i = 0 ; i < ed->GetNumberOfSteps(currentCalID) ; ++i)
        if ( validset.find(i) == validset.end())
            markLine(i);
    UpdateMe();
}

void EditorCanvas::markInterval(Int_t y)
{
    TH2* h = dynamic_cast<TH2*>(fClickSelected);
    if (h){
        auto index = floor(AbsPixeltoY(y));

        if (!intervalStartSet)
        {
            indexInterVal.Start() = index;
            fillLine(index);
            intervalStartSet = true;
        }
        else
        {
            indexInterVal.Stop() = index;
            indexInterVal.MakeSane();
            updateCalHist();
            for (auto inInt = indexInterVal.Start(); inInt <= indexInterVal.Stop(); ++inInt)
                fillLine(inInt);
            cout << "   Calibration steps marked:" << endl
                 << "     [" << indexInterVal.Start() << ", " << indexInterVal.Stop() <<"]" << endl
                 << endl;
        }
    }
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
        cout << "   Calibration steps marked:" << endl;
        for (const auto& in: indexMemory)
            cout << in << endl;
        cout << endl;
    }
}

void EditorCanvas::HandleInput(EEventType button, Int_t x, Int_t y)
{
    TCanvas::HandleInput(button,x,y);

    //switch (state)
    //{
    //case state_t::expand:
    //case state_t::remove:
    if (button == kButton1Down)
        markLine(y);
//    break;
//    case state_t::cut:
//        if (button == kButton1Down)
//            markInterval(y);
//        break;
//    default:
//        break;
//    }


}
