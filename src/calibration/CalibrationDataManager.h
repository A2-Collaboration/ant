#pragma once

//ant
#include "tree/TCalibrationData.h"
#include "tree/TDataRecord.h"
#include "base/std_ext.h"
#include "base/interval.h"

//std
#include <map>
#include <vector>
#include <list>
#include <string>
#include <memory>


namespace ant
{
namespace calibration
{

using DataMap_t = std::map<std::string,std::vector<TCalibrationData>>;

class DataBase
{
private:
    const std::string cm_treename_prefix;
    const std::string cm_branchname;

public:
    void ReadData(const std::string& filename);
    void WriteData(const std::string& filename);
    DataMap_t DataMap;
    DataBase():
        cm_treename_prefix("calibration-"),
        cm_branchname("cdata")
    {}
};

class DataAccess
{
public:
   /**
    *  \brief GetData Query the calibration database for specific TID
    *  \param calibrationID Calibration ID
    *  \param eventID event ID
    *  \param cdata   Reference to a TCalibrationData, data will be writter here
    *  \return true if valid data was found
    */
    virtual void Add(const TCalibrationData& data) = 0;

    virtual bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) = 0;

    virtual const std::list<TID> GetChangePoints(const std::string& calibrationID) = 0;

    virtual bool GetLastEntry(const std::string& calibrationID, TCalibrationData& cdata) = 0;
};


class DataManager: public DataAccess
{

private:

    std::string dataFileName;
    bool changedDataBase;
    std::unique_ptr<DataBase> dataBase;

    void lazyInit();

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

    /**
     * @brief getDepth returns the distance in steps to the latest calibration iteration
     * @param tid      event id
     * @param calibrationID  calibration id
     * @return depth
     */
    std::uint32_t getDepth(const TID& tid, const std::string& calibrationID) const;

public:
    DataManager(const std::string& DataFileName):
        dataFileName(DataFileName),
        changedDataBase(false)
    {}


    ~DataManager()
    {
        if (changedDataBase)
            dataBase->WriteData(dataFileName);
    }


    void Add(const TCalibrationData& data) override
    {
        lazyInit();
        dataBase->DataMap[data.CalibrationID].push_back(data);
        changedDataBase = true;
    }

    bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) override;

    const std::list<TID> GetChangePoints(const std::string& calibrationID) override;

    std::uint32_t GetNumberOfCalibrations()
    {
        lazyInit();
        return dataBase->DataMap.size();
    }

    std::uint32_t GetNumberOfDataPoints(const std::string& calibrationID) ;

    bool GetIDRange(const std::string& calibrationID, interval<TID>& IDinterval) ;

    bool GetLastEntry(const std::string& calibrationID, TCalibrationData& cdata) override;

};


} //namespace calibration
} //namespace ant
