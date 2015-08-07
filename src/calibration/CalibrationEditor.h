#pragma once
#include "analysis/plot/root_draw.h"
#include "CalibrationDataManager.h"

#include <memory>
#include <string>

namespace ant
{

namespace calibration
{

class Editor
{
private:
    DataBase dman;
    ant::canvas theCanvas;


    bool getIDRange(const std::string& calibrationID, interval<TID>& IDinterval) const;

public:
    Editor(): dman(){}

    void AddFromFile(const std::string& fileName) {dman.ReadData(fileName); }
    void SaveToFile(const std::string& fileName) const {dman.WriteData(fileName);}


    void ShowHistory(const std::string& calibrationID) const;


};

}
}
