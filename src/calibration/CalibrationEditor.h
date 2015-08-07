#pragma once
#include "analysis/plot/root_draw.h"
#include "CalibrationDataManager.h"

#include <memory>
#include <string>
#include <vector>

namespace ant
{

namespace calibration
{

class Editor
{
private:
    DataBase dman;

    bool getIDRange(const std::string& calibrationID, interval<TID>& IDinterval) const;
    std::vector<ant::TCalibrationData> getValidData(const std::string& calibrationID) const;

public:
    Editor(): dman(){}

    void AddFromFile(const std::string& fileName) {dman.ReadData(fileName); }
    void SaveToFile(const std::string& fileName) const {dman.WriteData(fileName);}

    void Add(const TCalibrationData& cdata) {dman.DataMap[cdata.CalibrationID].push_back(cdata);}

    void ListCalibrations() const;

    void ShowHistory(const std::string& calibrationID) const;

    /**
     * @brief ShowValid shows only calibration steps which are used
     * @param calibrationID
     */
    void ShowValid(const std::string&  calibrationID) const;



};

}
}
