#pragma once


#include "tree/TDataRecord.h"
#include "tree/TCalibrationData.h"
#include "base/interval.h"


#include <set>
#include <list>
#include <stdexcept>

namespace ant {

class WrapTFileOutput;

namespace calibration {


class DataBase
{

protected:

    std::string calibrationDataFolder;

    struct Range_t : interval<TID> {
        std::string FolderPath;
        bool operator<(const Range_t& other) const {
            return Start() < other.Start();
        }
        using interval<TID>::interval;
        Range_t(const TID& tid) : interval<TID>(tid, TID()) {}
        Range_t(const interval<TID> tidRange, const std::string& folderPath) :
            interval<TID>(tidRange),
            FolderPath(folderPath)
        {}
    };

    std::set<Range_t> getDataRanges(const std::string& calibrationID) const;

    bool loadFile(const std::string& filename, TCalibrationData& cdata) const;
    bool writeFile(const std::string& folder, const TCalibrationData& cdata) const;

    void handleStrictRange(const TCalibrationData& cdata) const;
    void handleRightOpen(const TCalibrationData& cdata) const;

    std::string makeTIDString(const TID& time) const;
    interval<TID> parseTIDRange(const std::string& tidRangeStr) const;
public:


    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

    std::list<std::string> GetCalibrationIDs() const;
    size_t GetNumberOfCalibrationData(const std::string& calibrationID) const;


    bool GetItem(const std::string& calibrationID,
                 const TID& currentPoint,
                 TCalibrationData& theData,
                 TID& nextChangePoint) const;

    enum class mode_t
    {
        AsDefault,
        RightOpen,
        StrictRange
    };
    void AddItem(const TCalibrationData& cdata, mode_t mode);


    DataBase(const std::string calibrationDataFolder_);
};

}//calibration
}//ant
