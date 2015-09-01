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
public:
    using DataMap_t = std::map<std::string,std::vector<ant::TCalibrationData>>;
private:
    const std::string cm_treename_prefix;
    const std::string cm_branchname;

    /**
     * @brief getDepth returns the distance in steps to the latest calibration iteration
     * @param tid      event id
     * @param calibrationID  calibration id
     * @return depth
     */
    std::uint32_t getDepth(const TID& tid, const std::string& calibrationID) const;

    /**
     * @brief isValid tests if the given id is a valid changepoint ( no newer data exists )
     * @param tid      queried event id
     * @param calibrationID  calibration id
     * @param depth    depth(distance from last calibration iteration) of given data point
     * @return valid or not
     */
    bool isValid(const TID& tid, const std::string& calibrationID, const std::uint32_t& depth) const
    {
        return (depth <= getDepth(tid,calibrationID));
    }

    void WriteToTree(WrapTFileOutput& file,
                     const DataMap_t::value_type& calibration) const;
    bool changed = false;
    DataMap_t dataMap;
public:
    /**
     * @brief ReadFromFile read all TCalibrationData trees from given filename
     * @param filename given ROOT file
     * @return true if successful
     */
    bool ReadFromFile(const std::string& filename);
    /**
     * @brief WriteToFile writes all data into given filename
     * @param filename
     */
    void WriteToFile(const std::string& filename) const;

    /**
     * @brief ReadFromFolder open all *.root and read TCalibrationData from them
     * @param folder
     * @return true if successful
     */
    bool ReadFromFolder(const std::string& folder);
    /**
     * @brief WriteToFolder write data into separate files, given by TCalibrationData::CalibrationID
     * @param folder where files will be created
     */
    void WriteToFolder(const std::string& folder) const;


    bool Has(const std::string& calibrationID) const { return dataMap.count(calibrationID) != 0; }
    std::uint32_t GetNumberOfDataPoints(const std::string& calibrationID) const;
    const std::list<TID> GetChangePoints(const std::string& calibrationID) const;

    const DataMap_t& DataMap() const { return dataMap; }
    DataMap_t::mapped_type& ModifyItem(const std::string& calibrationID) { return dataMap.at(calibrationID); }
    void AddCalibrationData(const TCalibrationData& cdata) { return dataMap[cdata.CalibrationID].push_back(cdata); }


    DataBase();
};

}//calibration
}//ant
