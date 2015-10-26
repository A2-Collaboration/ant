#pragma once

//std
#include <list>
#include <string>
#include <memory>


namespace ant
{

class TCalibrationData;
class TID;

namespace calibration
{

class DataBase;

class DataAccess
{
public:
    /**
     * @brief Add the given calibration data to the database
     * @param cdata
     */
    virtual void Add(const TCalibrationData& cdata) = 0;

    /**
    *  \brief GetData Query the calibration database for specific TID
    *  \param calibrationID Calibration ID
    *  \param eventID event ID
    *  \param cdata   Reference to a TCalibrationData, data will be writter here
    *  \return true if valid data was found
    */
    virtual bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) = 0;

    /**
     * @brief GetChangePoints obtains the IDs where the data should be changed
     * @param calibrationID the calibration identifier
     * @return list of change points (possibly unsorted)
     */
    virtual const std::list<TID> GetChangePoints(const std::string& calibrationID) = 0;
};


class DataManager: public DataAccess
{

private:

    std::string calibrationDataFolder;
    std::unique_ptr<DataBase> dataBase;
    bool extendable;

    void Init();



public:
    DataManager(const std::string& calibrationDataFolder_);


    ~DataManager();


    void Add(const TCalibrationData& cdata) override;

    bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) override;

    const std::list<TID> GetChangePoints(const std::string& calibrationID) override;

    std::uint32_t GetNumberOfCalibrations();

    std::uint32_t GetNumberOfDataPoints(const std::string& calibrationID);

    void SetExtendable() { extendable = true; }

};


} //namespace calibration
} //namespace ant
