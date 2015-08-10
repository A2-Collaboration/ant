#include "ACECanvas.h"
#include "calibration/CalibrationEditor.h"
#include "TH2D.h"
#include "TH2.h"


#include <cmath>

#include <iostream>

using namespace ant;
using namespace std;


void ACECanvas::update_modified()
{
    Modified();
    Update();
}

void ACECanvas::makeCalHist()
{
    calHist = new TH2D( (std_ext::formatter() << "hist-" << currentCalID).str().c_str(),
                        (std_ext::formatter() << "History for " << currentCalID).str().c_str(),
                        100, 0, 100,
                        ed.GetNumberOfSteps(currentCalID),
                        0, ed.GetNumberOfSteps(currentCalID)
                        );
    calHist->SetXTitle("TID [%]");
    calHist->SetYTitle("Calibration Step");
    calHist->Draw("col");
}

void ant::ACECanvas::loadFile(const std::string& fileName)
{
    //    if ( !ed->AddFromFile(fileName) )
    //        throw runtime_error(string("Could not file ")+filename);

    //DEBUG!!!!
    cout << "try load " << fileName << "  in future!!!" << endl;
    addSomeRandomData();
    currentCalID = "testID";
    // end DEBUG!!!!!!
    makeCalHist();
    updateCalHist();

    change_state(state_t::base);
}

void ant::ACECanvas::updateCalHist()
{
    calHist->Reset();
    for (const auto& ran: ed.GetAllRanges(currentCalID))
        for (int i = floor(ran.second.Start()); i < ceil(ran.second.Stop()) ; ++i)
            calHist->SetBinContent(i+1,ran.first+1,2);

    for (const auto& ran: ed.GetAllValidRanges(currentCalID))
        for (int i = floor(ran.second.Start()); i < ceil(ran.second.Stop()) ; ++i)
            calHist->SetBinContent(i+1,ran.first+1,1);

    update_modified();
}

//void ant::TACECanvas::SaveToFile(const std::string& fileName)
//{
//    ed->SaveToFile(fileName);
//}

void ACECanvas::change_state(ACECanvas::state_t newstate)
{
    switch (newstate) {
    cout << endl;
    case state_t::base:
        cout << "Base Mode:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  r   remove Calibration steps" << endl
             << "   *  c   cut out intervals of Calibration steps" << endl;
        break;
    case state_t::remove:
        cout << "Remove Mode:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  a   apply changes" << endl
             << "   *  c   cancel" << endl;
        break;
    case state_t::cut:
        cout << "Cut interval Mode:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  a   apply changes" << endl
             << "   *  c   cancel" << endl;
        break;
    default:
        break;
    }
    cout << endl;
    state = newstate;
}

ACECanvas::ACECanvas(const string &FileName):
    TCanvas("Ant Calibration Editor"),
    state(state_t::init),
    fileName(FileName),
    currentCalID(),
    ed(),
    intervalStartSet(false)
{
    loadFile(fileName);
}

void ACECanvas::removeAllinStepMemory()
{
    cout << endl << "Removing:" << endl;
    for (const auto& in: indexMemory)
        if (ed.Remove(currentCalID,in))
            cout << "  " << in << endl;
    indexMemory.clear();
    updateCalHist();
}
void ACECanvas::removeInterValFromMemory()
{
    if (!intervalStartSet)
        return;
    cout << endl << "Removing:" << endl;
    if (ed.Remove(currentCalID,indexInterVal[0],indexInterVal[1]))
        cout << "  [" << indexInterVal[0] << ", " << indexInterVal[1] << "]" << endl;
    intervalStartSet = false;
    updateCalHist();
}

void ACECanvas::HandleKeypress(const char key)
{
    switch (state)
    {
    case state_t::base:
        switch (key)
        {
        case 'r':
            change_state(state_t::remove);
            break;
        case 'c':
            change_state(state_t::cut);
            break;
        default:
            break;
        }
        break;
    case state_t::remove:
        switch(key)
        {
        case 13:
        case 'a':
            removeAllinStepMemory();
            change_state(state_t::base);
            break;
        case 'c':
            indexMemory.clear();
            updateCalHist();
            change_state(state_t::base);
            break;
        default:
            break;
        }
        break;
    case state_t::cut:
        switch(key)
        {
        case 13:
        case 'a':
            removeInterValFromMemory();
            change_state(state_t::base);
            break;
        case 'c':
            intervalStartSet = false;
            updateCalHist();
            change_state(state_t::base);
            break;
        default:
            break;
        }
        break;
    default:
        break;
    }

}

void ACECanvas::fillLine(uint32_t lineNumber)
{
    for (Int_t i = 0; i < 100; ++i )
        calHist->Fill(i,lineNumber,3);
}

void ACECanvas::unFillLine(uint32_t lineNumber)
{
    for (Int_t i = 0; i < 100; ++i )
        calHist->Fill(i,lineNumber,-3);
}

void ACECanvas::markInterval(Int_t y)
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
            cout << "   Calibration steps marked for remove:" << endl
                 << "     [" << indexInterVal.Start() << ", " << indexInterVal.Stop() <<"]" << endl
                 << endl;
        }
    }
}

void ACECanvas::markLine(Int_t y)
{
    TH2* h = dynamic_cast<TH2*>(fClickSelected);
    if (h){
        auto step = floor(AbsPixeltoY(y));
        auto foundat = indexMemory.find(step);
        if ( foundat == indexMemory.end() )
        {
            indexMemory.emplace(step);
            fillLine(step);
//            for (Int_t i = minx; i < maxx; ++i )
//                calHist->Fill(i-1,step-1,3);
        }
        else
        {
            indexMemory.erase(foundat);
            unFillLine(step);
//            for (Int_t i = minx; i < maxx; ++i )
//                calHist->Fill(i-1,step-1,-3);
        }
        cout << "   Calibration steps marked for remove:" << endl;
        for (const auto& in: indexMemory)
            cout << in << endl;
        cout << endl;
    }
//    auto bin = h->GetYaxis()->FindBin(y);
}

void ACECanvas::HandleInput(EEventType button, Int_t x, Int_t y)
{
    TCanvas::HandleInput(button,x,y);

    switch (state)
    {
    case state_t::remove:
        if (button == kButton1Down)
            markLine(y);
        break;
    case state_t::cut:
        if (button == kButton1Down)
            markInterval(y);
        break;
    default:
        break;
    }


    if(button == kKeyPress)
        HandleKeypress(x);
}

//DEBUG
void ACECanvas::addSomeRandomData()
{
    unsigned time = 0;
    auto makedata = [&time] (unsigned first, unsigned last)
    {
        TCalibrationData tmp("Wolfes",
                             "testData",
                             time++,
                             "testID",
                             TID(0,first),TID(0,last));
        tmp.Data.emplace_back(0,1);
        tmp.Data.emplace_back(1,2);
        return tmp;
    };

    ed.Add(makedata(0,1000));
    ed.Add(makedata(0,100));
    ed.Add(makedata(450,460));
    ed.Add(makedata(100,200));
    ed.Add(makedata(12,200));
    ed.Add(makedata(420,900));
    for (int i = 50 ; i < 100; ++i)
        ed.Add(makedata(i*10,(i+1)*10));
    for (int i = 50 ; i < 100; ++i)
        ed.Add(makedata(i*10,(i+1)*10));
    ed.Add(makedata(800,900));
}
