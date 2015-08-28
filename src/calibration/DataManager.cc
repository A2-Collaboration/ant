#include "DataManager.h"

//ant
#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/memory.h"
#include "base/interval.h"

//ROOT
#include "TTree.h"
#include "TError.h"
#include "TDirectory.h"

//std
#include <algorithm>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::calibration;


void DataManager::lazyInit()
{
    if ( dataBase == nullptr )
    {
        dataBase = std_ext::make_unique<DataBase>();
        dataBase->ReadData(dataFileName);
    }
}

void DataManager::Add(const TCalibrationData& data)
{
    lazyInit();
    dataBase->DataMap[data.CalibrationID].push_back(data);
    changedDataBase = true;
}


bool DataManager::GetData(const string& calibrationID, const TID& eventID, TCalibrationData& cdata)
{
    lazyInit();
    //case one: calibration doesn't exist
    if ( dataBase->DataMap.count(calibrationID) == 0)
        return false;

    //case two: calibration exists
    auto& calibPairs = dataBase->DataMap.at(calibrationID);
    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
    {
        interval<TID> range(rit->FirstID,rit->LastID);
        if (range.Contains(eventID))
        {
            cdata = *rit;
            return true;
        }
    }

    //case three: TID not covered by calibration
    return false;
}


const list<TID> DataManager::GetChangePoints(const string& calibrationID)
{
    lazyInit();
    return dataBase->GetChangePoints(calibrationID);
}

uint32_t DataManager::GetNumberOfCalibrations()
{
    lazyInit();
    return dataBase->DataMap.size();
}

uint32_t DataManager::GetNumberOfDataPoints(const string& calibrationID)
{
    lazyInit();
    return dataBase->GetNumberOfDataPoints(calibrationID);
}



bool DataManager::GetLastEntry(const std::string& calibrationID, TCalibrationData& cdata)
{
    lazyInit();
    if (dataBase->DataMap.count(calibrationID) == 0)
        return false;

    auto& data = dataBase->DataMap.at(calibrationID);

    cdata = data.back();

    return true;
}



