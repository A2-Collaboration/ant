#pragma once

#include "Calibration.h"

//std
#include <list>
#include <string>
#include <memory>

namespace ant
{

struct TCalibrationData;
struct TID;

namespace calibration
{

class DataAccess
{
public:
    /**
     * @brief Add the given calibration data to the database
     * @param cdata
     */
    virtual void Add(const TCalibrationData& cdata, Calibration::AddMode_t addMode) = 0;

    /**
    *  \brief GetData Query the calibration database for specific TID
    *  \param calibrationID Calibration ID
    *  \param eventID event ID
    *  \param cdata   Reference to a TCalibrationData, data will be writter here
    *  \return true if valid data was found
    */
    virtual bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) = 0;

};

class DataBase;

class DataManager : public DataAccess
{

private:

    const std::string calibrationDataFolder;
    std::unique_ptr<DataBase> dataBase;

    void Init();

    bool override_as_default = false;

public:
    DataManager(const std::string& calibrationDataFolder_);
    virtual ~DataManager();

    void Add(const TCalibrationData& cdata, Calibration::AddMode_t addMode) override;

    bool GetData(const std::string& calibrationID,
                 const TID& eventID,
                 TCalibrationData& cdata) override;
    bool GetData(const std::string& calibrationID,
                 const TID& eventID,
                 TCalibrationData& cdata,
                 TID& nextChangePoint);

    // the following methods are only useful for test cases
    std::list<std::string> GetCalibrationIDs();
    std::size_t GetNumberOfCalibrationIDs();
    std::size_t GetNumberOfCalibrationData(const std::string& calibrationID);

    bool GetOverrideToDefault() const;
    void SetOverrideToDefault(bool v);

    std::string GetCalibrationDataFolder() const;

};


} //namespace calibration
} //namespace ant
