#include "ACECanvas.h"
#include "calibration/CalibrationEditor.h"
#include "TH2D.h"

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

    change_state(state_t::base);
}

//void ant::TACECanvas::SaveToFile(const std::string& fileName)
//{
//    ed->SaveToFile(fileName);
//}

void ACECanvas::change_state(ACECanvas::state_t newstate)
{
    switch (newstate) {
    case state_t::base:
//        updateHist();
        break;
    default:
        break;
    }

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

void ACECanvas::HandleKeypress(const char key)
{
    switch (key) {
    default:
        break;
    }
}
void ACECanvas::HandleInput(EEventType button, Int_t x, Int_t y)
{
    TCanvas::HandleInput(button,x,y);

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
