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


void DataManager::Init()
{
    if ( dataBase == nullptr )
    {
        dataBase = std_ext::make_unique<DataBase>();
        dataBase->ReadFromFolder(calibrationDataFolder);
    }
}

void DataManager::Add(const TCalibrationData& cdata)
{
    Init();
    dataBase->AddItem(cdata);
    changedDataBase = true;
}


bool DataManager::GetData(const string& calibrationID, const TID& eventID, TCalibrationData& cdata)
{
    Init();
    //case one: calibration doesn't exist
    if ( dataBase->DataMap().count(calibrationID) == 0)
        return false;

    //case two: calibration exists
    const auto& calibPairs = dataBase->GetItems(calibrationID);
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
    Init();
    return dataBase->GetChangePoints(calibrationID);
}

uint32_t DataManager::GetNumberOfCalibrations()
{
    Init();
    return dataBase->DataMap().size();
}

uint32_t DataManager::GetNumberOfDataPoints(const string& calibrationID)
{
    Init();
    return dataBase->GetNumberOfDataPoints(calibrationID);
}



bool DataManager::GetLastEntry(const std::string& calibrationID, TCalibrationData& cdata)
{
    Init();
    if (dataBase->DataMap().count(calibrationID) == 0)
        return false;

    const auto& data = dataBase->GetItems(calibrationID);

    cdata = data.back();

    return true;
}



