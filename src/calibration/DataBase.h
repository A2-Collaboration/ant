#pragma once


#include "tree/TCalibrationData.h"
#include "base/interval.h"
#include "Calibration.h"

#include <list>
#include <stdexcept>
#include <memory>

namespace ant {

class WrapTFileOutput;

namespace calibration {


class DataBase
{
public:

    DataBase(const std::string& calibrationDataFolder);

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

    bool GetItem(const std::string& calibrationID,
                 const TID& currentPoint,
                 TCalibrationData& theData,
                 TID& nextChangePoint) const;

    /**
     * @brief Load TObject from current calibration ROOT file
     * @param calibrationID   ID of the calibration to access
     * @param currentPoint    TID of current event
     * @param objname         ROOT name of the TObejct to look for
     * @param nextChangePoint
     * @return A Clone of the TObject in the file, nullptr if not found.
     */
    TObject* GetTObject(const std::string& calibrationID,
                 const TID& currentPoint,
                 const std::string& objname,
                 TID& nextChangePoint) const;

    void AddItem(const TCalibrationData& cdata, Calibration::AddMode_t mode);

    std::list<std::string> GetCalibrationIDs() const;
    size_t GetNumberOfCalibrationData(const std::string& calibrationID) const;

    struct OnDiskLayout {

        /**
         * @brief EnableCaching if true, the OnDiskLayout does not scan the folder structure
         * every time GetDataRanges is called. Note that this should only be enabled globally if
         * only read accesses are executed. The cache prevents changes to be seen made by new items!
         */
        static bool EnableCaching;

        const std::string CalibrationDataFolder;

        OnDiskLayout(const std::string& calibrationDataFolder);

        enum class Type_t {
            DataDefault, DataRanges, MC
        };

        std::string GetFolder(const std::string& calibrationID, Type_t type) const;
        std::string GetCurrentFile(const std::string& calibrationID, Type_t type) const;
        std::string GetRangeFolder(const std::string& calibrationID, const interval<TID>& range) const;
        std::string RemoveCalibrationDataFolder(const std::string& path) const;

        struct Range_t : interval<TID> {
            std::string FolderPath;
            using interval<TID>::interval;
            bool operator<(const Range_t& other) const {
                return Start() < other.Start();
            }
            Range_t(const interval<TID> tidRange, const std::string& folderPath) :
                interval<TID>(tidRange),
                FolderPath(folderPath)
            {}
        };
        std::string GetCurrentFile(const Range_t& range) const;
        using DataRanges_t = std::list<Range_t>;
        DataRanges_t GetDataRanges(const std::string& calibrationID) const;

    protected:
        std::string makeTIDString(const TID& tid) const;
        interval<TID> parseTIDRange(const std::string& tidRangeStr) const;
        mutable std::map<std::string, DataRanges_t> cached_ranges;
    };



protected:
    OnDiskLayout Layout;

    bool loadFile(const std::string& filename, TCalibrationData& cdata) const;
    TObject* loadObjectFromFile(const std::string& filename, const std::string& objectname) const;

    bool writeToFolder(const std::string& folder, const TCalibrationData& cdata) const;

    void handleStrictRange(const TCalibrationData& cdata) const;
    void handleRightOpen(const TCalibrationData& cdata) const;
};

}//calibration
}//ant
