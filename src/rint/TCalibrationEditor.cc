#include "TCalibrationEditor.h"
#include "calibration/CalibrationEditor.h"

#include <iostream>

using namespace ant;
using namespace std;

ant::TCalibrationEditor::TCalibrationEditor()
{
    ed = new calibration::Editor();
}


void ant::TCalibrationEditor::AddFromFile(const std::string& fileName)
{
    ed->AddFromFile(fileName);
}

void ant::TCalibrationEditor::SaveToFile(const std::string& fileName)
{
    ed->SaveToFile(fileName);
}

void TCalibrationEditor::Remove(const std::string& calibrationID, const uint32_t& index)
{
    if (ed->Remove(calibrationID,index))
        cout << "Succesfully removed iteration " << index << " in calibration " << calibrationID << "." << endl;
    else
        cout << "Calibration " << calibrationID << "or stepnumber don't exist." << endl;
}

void TCalibrationEditor::Remove(const std::string& calibrationID, const uint32_t& index1, const uint32_t& index2)
{
    if (ed->Remove(calibrationID,index1,index2))
        cout << "Succesfully removed iterations [" << index1 << ", " << index2 << "]"
             << " in calibration " << calibrationID << "." << endl;
    else
        cout << "Calibration " << calibrationID << " or stepnumbers don't exist." << endl;
}

void TCalibrationEditor::ReduceToValid(const string &calibrationID)
{
    if (ed->ReduceToValid(calibrationID))
        cout << "Succesfully reduced calibration " << calibrationID << "." << endl;
    else
        cout << "Unable to reduce caibration " << calibrationID << "." << endl;
}

void ant::TCalibrationEditor::ShowHistory(const std::string& calibrationID) const
{
    ed->ShowHistory(calibrationID);
}

void TCalibrationEditor::ShowValid(const std::string& calibrationID) const
{
    ed->ShowValid(calibrationID);
}

ant::TCalibrationEditor::~TCalibrationEditor()
{
    delete ed;
}

void TCalibrationEditor::AddSomeRandomData()
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

    ed->Add(makedata(0,1000));
    ed->Add(makedata(0,100));
    ed->Add(makedata(450,460));
    ed->Add(makedata(100,200));
    ed->Add(makedata(12,200));
    ed->Add(makedata(420,900));
    for (int i = 50 ; i < 100; ++i)
        ed->Add(makedata(i*10,(i+1)*10));
    for (int i = 50 ; i < 100; ++i)
        ed->Add(makedata(i*10,(i+1)*10));
    ed->Add(makedata(800,900));
}

void TCalibrationEditor::ListCalibrations() const
{
    ed->ListCalibrations();
}
