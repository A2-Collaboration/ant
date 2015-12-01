#pragma once

#include <map>
#include <vector>
#include <list>

#include "tree/TDataRecord.h"
#include "tree/TCalibrationData.h"


namespace ant {

class WrapTFileOutput;

namespace calibration {


class DataBase
{

protected:

    std::string calibrationDataFolder;

public:
    enum class mode_t
    {
        AsDefault,
        RightOpen,
        StrictRange
    };

    bool Has(const std::string& calibrationID) const;
    std::list<std::string> GetKeys() const;
    std::uint32_t GetNumberOfDataPoints(const std::string& calibrationID) const;


    bool GetItem   (const std::string& calibrationID, const TID& currentPoint,TCalibrationData& theData, TID& nextChangePoint) const;

    void AddItem(const TCalibrationData& cdata, mode_t mode);


    DataBase(const std::string calibrationDataFolder_);
};

}//calibration
}//ant
