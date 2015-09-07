#include "ACECanvas.h"

#include "calibration/Editor.h"
#include "base/std_ext/string.h"

#include "TH2D.h"
#include "TH2.h"

#include <cmath>

#include <iostream>

using namespace ant;
using namespace std;
using namespace ant::calibration::gui;


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
    calHist->SetStats(false);
}

void ACECanvas::loadFile(const std::string& fileName)
{
    if ( !ed.AddFromFile(fileName) )
        throw runtime_error(string("Could not file ")+fileName);

    if (ed.GetListOfCalibrations().size() == 0 )
        throw runtime_error(string("No calibration found in file ")+fileName);

    cout << endl
         << "  Sucessfully opened file."<< endl;

//    loadCalibration();

//    if (currentCalID.empty())
//        throw runtime_error("No calibration loaded, exiting... ");

    change_state(state_t::base);
}

void ACECanvas::loadCalibration()
{
    changeCalibrationID();

    makeCalHist();
    updateCalHist();
}

void ACECanvas::changeCalibrationID()
{
    cout << endl
         << "  Available calibrationIDs:" << endl;
    for (const auto& cID: ed.GetListOfCalibrations())
    {
        cout << "    *  " << cID << endl;
    }
    cout << endl
         << "Enter the calibration you want to edit," << endl
         << "'.q' to abort." << endl
         << "   -->  ";

    bool found = false;
    string calID;
    while (!found)
    {
        cin >> calID;
        if ( calID.compare(".q") == 0 )
            return;
        found = ed.Has(calID);
    }
    currentCalID = calID;
}

void ACECanvas::updateCalHist()
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

void ACECanvas::saveToFile(const std::string& fileName)
{
    ed.SaveToFile(fileName);
    cout << endl << "Saved to File " << fileName << "." << endl;
}

void ACECanvas::resetCalibrationTo(const string& calID)
{
    currentCalID = calID;
    intervalStartSet = false;
    indexMemory.clear();
    makeCalHist();
    updateCalHist();
}

void ACECanvas::change_state(ACECanvas::state_t newstate)
{
    cout << endl
         << endl
         << "=====  Ant Calibration Editor  ========================" << endl
         << endl;
    switch (newstate) {
    case state_t::init:
        cout << "Init: Choose Calibration" << endl;
        resetCalibrationTo("testID");
//        openQuery();
        change_state(state_t::base);
        break;
    case state_t::base:
        cout << "Base Mode:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  e   Expand Calibration steps" << endl
             << "   *  c   Cut out intervals of Calibration steps" << endl
             << "   *  l   open a load new CalibrationID dialog" << endl
             << "   *  r   Remove Calibration steps" << endl
             << "   *  s   Save Your changes" << endl
             << "   *  v   reduce to Valid Calibration steps" << endl;
        break;
    case state_t::save:
        cout << "Save all changes:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  o   Override existing file with changes" << endl
             << "   *  a   save As a new file" << endl
             << "   *  c   Cancel" << endl;
        break;
    case state_t::expand:
        cout << "Expand Mode:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  a   Apply changes" << endl
             << "   *  c   Cancel" << endl;
        break;
    case state_t::remove:
        cout << "Remove Mode:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  a   Apply changes" << endl
             << "   *  c   Cancel" << endl;
        break;
    case state_t::cut:
        cout << "Cut interval Mode:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  a   Apply changes" << endl
             << "   *  c   Cancel" << endl;
        break;
    case state_t::reduceToValid:
        markUnValid();
        cout << "Reduce Mode:" << endl
             << "  Keys in Canvas:" << endl
             << "   *  a   Apply changes" << endl
             << "   *  c   Cancel" << endl;
        break;
    default:
        break;
    }
    cout << endl
         << "=======================================================" << endl
         << endl;
    state = newstate;
}

ACECanvas::ACECanvas(const string &FileName):
    TCanvas("Ant Calibration Editor"),
    state(state_t::init),
    fileName(FileName),
    currentCalID(),
    ed(),
    selector(nullptr),
    intervalStartSet(false)
{
    loadFile(fileName);
}

void ACECanvas::expandAllinIndexMemory()
{
    cout << endl << "Expanding:" << endl;
    for (const auto& in: indexMemory)
        if (ed.ExpandToMax(currentCalID,in))
            cout << "  " << in << endl;
    indexMemory.clear();
    updateCalHist();
}
void ACECanvas::removeAllinIndexMemory()
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
void ACECanvas::openQuery()
{
    new ListQuery(this,"title","text",ed.GetListOfCalibrations());
}

void ACECanvas::HandleKeypress(const char key)
{
    switch (state)
    {
    case state_t::base:
        switch (key)
        {
        case 'e':
            change_state(state_t::expand);
            break;
        case 'l':
            openQuery();
            break;
        case 's':
            change_state(state_t::save);
            break;
        case 'r':
            change_state(state_t::remove);
            break;
        case 'c':
            change_state(state_t::cut);
            break;
        case 'v':
            change_state(state_t::reduceToValid);
            break;
        default:
            break;
        }
        break;
    case state_t::save:
        switch(key)
        {
        case 'o':
            saveToFile(fileName);
            change_state(state_t::base);
            break;
//        case 'a':
//            indexMemory.clear();
//            updateCalHist();
//            change_state(state_t::base);
//            break;
        case 'c':
            change_state(state_t::base);
            break;
        default:
            break;
        }
        break;
    case state_t::expand:
        switch(key)
        {
        case 13:
        case 'a':
            expandAllinIndexMemory();
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
    case state_t::remove:
        switch(key)
        {
        case 13:
        case 'a':
            removeAllinIndexMemory();
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
    case state_t::reduceToValid:
        switch(key)
        {
        case 13:
        case 'a':
            ed.ReduceToValid(currentCalID);
            updateCalHist();
            change_state(state_t::base);
            break;
        case 'c':
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

void ACECanvas::markUnValid()
{
    auto valids = ed.GetAllValidRanges(currentCalID);
    set<uint32_t> validset;
    for (const auto v: valids)
        validset.emplace(v.first);

    for (uint32_t i = 0 ; i < ed.GetNumberOfSteps(currentCalID) ; ++i)
        if ( validset.find(i) == validset.end())
            fillLine(i);
    update_modified();
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
            cout << "   Calibration steps marked:" << endl
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

void ACECanvas::HandleInput(EEventType button, Int_t x, Int_t y)
{
    TCanvas::HandleInput(button,x,y);

    switch (state)
    {
    case state_t::expand:
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

void ACECanvas::captureReturnValue(Query* dialog, std::vector<string>& returnValue)
{
    if (dialog)
        cout << endl
             << "Changing Runset to " << returnValue.at(0) << "." << endl;
    resetCalibrationTo(returnValue.at(0));
}

//DEBUG
void ACECanvas::addSomeRandomData()
{
    unsigned time = 0;
    auto makedata = [&time] (unsigned first, unsigned last)
    {
        TCalibrationData tmp("Wolfes",
                             time++,
                             "testID",
                             TID(first),TID(last));
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


TextPad::TextPad(double xmargin, double ymargin):
    xMargin(xmargin), yMargin(ymargin)
{
    paveText = new TPaveText(xMargin,yMargin,1-xMargin,1-yMargin,"nb");
}

TextPad::TextPad(std::list<string>& lines, double xmargin, double ymargin):
    xMargin(xmargin),
    yMargin(ymargin)
{
    paveText = new TPaveText(xMargin,yMargin,1-xMargin,1-yMargin,"nb");
    for (const auto& line: lines)
        Add(line);
}

void TextPad::Draw(const string& option) const
{
    paveText->Draw(option.c_str());
}
