#pragma once


#include "tree/TDataRecord.h"
#include "tree/TCalibrationData.h"
#include "base/interval.h"


#include <set>
#include <list>

namespace ant {

class WrapTFileOutput;

namespace calibration {


class DataBase
{

protected:

    std::string calibrationDataFolder;

    struct Range_t {
        interval<TID> Range;
        std::string Filename;
        bool operator<(const Range_t& other) const {
            return Range.Start() < other.Range.Start();
        }
        Range_t(const TID& tid) : Range(tid, TID()) {}
    };

    std::set<Range_t> getRanges(const std::string& calibrationID) const;

    bool loadFile(const std::string& filename, TCalibrationData& cdata) const;

public:
    enum class mode_t
    {
        AsDefault,
        RightOpen,
        StrictRange
    };

    //bool Has(const std::string& calibrationID) const;
    std::list<std::string> GetKeys() const;
    std::uint32_t GetNumberOfDataItems(const std::string& calibrationID) const;


    bool GetItem(const std::string& calibrationID,
                 const TID& currentPoint,
                 TCalibrationData& theData,
                 TID& nextChangePoint) const;

    void AddItem(const TCalibrationData& cdata, mode_t mode);


    DataBase(const std::string calibrationDataFolder_);
};

}//calibration
}//ant
