#include "ACECanvas.h"
#include "calibration/CalibrationEditor.h"
#include "TH2D.h"
#include "TH2.h"


#include <cmath>

#include <iostream>

using namespace ant;
using namespace std;


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
            calHist->SetBinContent(i,ran.first+1,2);

    for (const auto& ran: ed.GetAllValidRanges(currentCalID))
        for (int i = floor(ran.second.Start()); i < ceil(ran.second.Stop()) ; ++i)
            calHist->SetBinContent(i,ran.first+1,1);

    this->cd();
    calHist->Draw("col");
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
             << "   *  r         remove Calibration steps" << endl
             << "   *  c   cut out intervals of Calibration steps" << endl;
        break;
    case state_t::remove:
        cout << "Remove Mode:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  <return>  remove selection" << endl
             << "   *  <esc>     abort" << endl;
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
    ed()
{
    loadFile(fileName);
}

void ACECanvas::removeAllinStepMemory()
{
    cout << endl << "Removing:" << endl;
    for (const auto& in: stepMemory)
        if (ed.Remove(currentCalID,in-1))
            cout << "  " << in-1 << endl;
    stepMemory.clear();
    updateCalHist();
    this->Update();
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

        default:
            break;
        }
    case state_t::remove:
    case state_t::cut:
        switch(key)
        {
        case 13:
            removeAllinStepMemory();
            change_state(state_t::base);
            break;
        case 27:
            stepMemory.clear();
            updateCalHist();
            change_state(state_t::base);
        default:
            cout << "unknown key in remove" << endl;
            break;
        }
    default:
        break;
    }

}

void ACECanvas::markLine(Int_t y)
{
    TH2* h = dynamic_cast<TH2*>(fClickSelected);
    if (h){
        auto step = floor(AbsPixeltoY(y))+1;

        Int_t minx = floor(ed.GetRange(currentCalID,step-1).second.Start());
        Int_t maxx = ceil(ed.GetRange(currentCalID,step-1).second.Stop());

        auto foundat = stepMemory.find(step);
        if ( foundat == stepMemory.end() )
        {
            stepMemory.emplace(step);
            for (Int_t i = minx; i < maxx; ++i )
                calHist->Fill(i-1,step-1,3);
        }
        else
        {
            stepMemory.erase(foundat);
            for (Int_t i = minx; i < maxx; ++i )
                calHist->Fill(i-1,step-1,-3);
        }
        cout << "   Calibration steps marked for remove:" << endl;
        for (const auto& in: stepMemory)
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
        {
            markLine(y);
        }
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
