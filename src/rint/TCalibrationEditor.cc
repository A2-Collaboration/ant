#include "TCalibrationEditor.h"
#include "calibration/CalibrationEditor.h"

using namespace ant;

ant::TCalibrationEditor::TCalibrationEditor()
{
    ed = new calibration::Editor();
}

void ant::TCalibrationEditor::AddFromFile(const std::string& fileName) {ed->AddFromFile(fileName);}

void ant::TCalibrationEditor::SaveToFile(const std::string& fileName) {ed->SaveToFile(fileName);}

void ant::TCalibrationEditor::ShowHistory(const std::string& calibrationID) const
{
    ed->ShowHistory(calibrationID);
}

ant::TCalibrationEditor::~TCalibrationEditor()
{
   delete ed;
}
