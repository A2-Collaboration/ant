#pragma once

//ant
#include "tree/TCalibrationData.h"
#include "tree/TDataRecord.h"
#include "base/std_ext.h"
#include "base/interval.h"
#include "DataBase.h"

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


    void Add(const TCalibrationData& data) override;

    bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) override;

    const std::list<TID> GetChangePoints(const std::string& calibrationID) override;

    std::uint32_t GetNumberOfCalibrations();

    std::uint32_t GetNumberOfDataPoints(const std::string& calibrationID);


    bool GetLastEntry(const std::string& calibrationID, TCalibrationData& cdata) override;

};


} //namespace calibration
} //namespace ant
