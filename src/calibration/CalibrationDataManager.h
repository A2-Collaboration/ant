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

class CalibrationDataManager
{
    using dataMap = std::map<std::string,std::vector<TCalibrationData>>;
private:
    const std::string cm_treename_prefix;
    const std::string cm_branchname;
    std::string dataFileName;
    bool changedDataBase;

    std::unique_ptr<dataMap> dataBase;

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

    /**
     * @brief finish takes care of rewriting the data to the tree
     */
    void writeDataBase() ;
    void readDataBase();

public:
    CalibrationDataManager(const std::string& DataFileName);

    ~CalibrationDataManager()
    {
        if (changedDataBase)
            writeDataBase();
    }


    void Add(const TCalibrationData& data)
    {
        lazyInit();
        dataBase->operator[](data.CalibrationID).push_back(data);
        changedDataBase = true;
    }

    /**
     *  \brief GetData Query the calibration database for specific TID
     *  \param calibrationID Calibration ID
     *  \param eventID event ID
     *  \param cdata   Reference to a TCalibrationData, data will be writter here
     *  \return true if valid data was found
     */
    bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) ;

    const std::list<TID> GetChangePoints(const std::string& calibrationID) ;

    std::uint32_t GetNumberOfCalibrations()
    {
        lazyInit();
        return dataBase->size();
    }

    std::uint32_t GetNumberOfDataPoints(const std::string& calibrationID) ;

    bool GetIDRange(const std::string& calibrationID, interval<TID>& IDinterval) ;

    bool GetLastEntry(const std::string& calibrationID, TCalibrationData& cdata) ;

};


} //namespace ant
