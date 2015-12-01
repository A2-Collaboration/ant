#pragma once

//std
#include <list>
#include <string>
#include <memory>

#include "calibration/DataBase.h"


namespace ant
{

class TCalibrationData;
class TID;

namespace calibration
{


class DataAccess
{
public:
    /**
     * @brief Add the given calibration data to the database
     * @param cdata
     */
    virtual void Add(const TCalibrationData& cdata, DataBase::mode_t addMode) = 0;

    /**
    *  \brief GetData Query the calibration database for specific TID
    *  \param calibrationID Calibration ID
    *  \param eventID event ID
    *  \param cdata   Reference to a TCalibrationData, data will be writter here
    *  \return true if valid data was found
    */
    virtual bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) = 0;

};


class DataManager : public DataAccess
{

private:

    std::string calibrationDataFolder;
    std::unique_ptr<DataBase> dataBase;

    void Init();



public:
    DataManager(const std::string& calibrationDataFolder_);


    ~DataManager();


    void Add(const TCalibrationData& cdata, DataBase::mode_t addMode) override;

    bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) override;
    bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata, TID& nextChangePoint);

    std::uint32_t GetNumberOfCalibrations();

    std::uint32_t GetNumberOfDataPoints(const std::string& calibrationID);

};


} //namespace calibration
} //namespace ant
