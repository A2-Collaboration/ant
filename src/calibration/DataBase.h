#pragma once

#include <map>
#include <vector>
#include <list>

#include "tree/TDataRecord.h"
#include "tree/TCalibrationData.h"


namespace ant
{
namespace calibration
{

using DataMap_t = std::map<std::string,std::vector<ant::TCalibrationData>>;

class DataBase
{
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

public:
    bool ReadData(const std::string& filename);
    void WriteData(const std::string& filename) const;

    bool Has(const std::string& calibrationID) const{ return (DataMap.count(calibrationID) != 0);}
    std::uint32_t GetNumberOfDataPoints(const std::string& calibrationID) const;
    const std::list<TID> GetChangePoints(const std::string& calibrationID) const;


    DataMap_t DataMap;
    DataBase():
        cm_treename_prefix("calibration-"),
        cm_branchname("cdata")
    {}
};

}//calibration
}//ant
