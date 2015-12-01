#include "DataManager.h"

//ant
#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/memory.h"
#include "base/interval.h"
#include "DataBase.h"
#include "tree/TCalibrationData.h"
#include "tree/TDataRecord.h"

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
    if(dataBase)
        return;
    dataBase = std_ext::make_unique<DataBase>();
    dataBase->ReadFromFolder(calibrationDataFolder);
}

DataManager::DataManager(const string& calibrationDataFolder_):
    calibrationDataFolder(calibrationDataFolder_)
{}

DataManager::~DataManager()
{
    if(dataBase)
        dataBase->WriteToFolder(calibrationDataFolder);
}

void DataManager::Add(const TCalibrationData& cdata, DataBase::mode_t addMode)
{
    Init();
    dataBase->AddItem(cdata, addMode);
    LOG(INFO) << "Added " << cdata;
}

bool DataManager::GetData(const string& calibrationID, const TID& eventID, TCalibrationData& cdata)
{
    TID dummy;
    return GetData(calibrationID, eventID, cdata, dummy);
}


bool DataManager::GetData(const string& calibrationID,
                          const TID& eventID, TCalibrationData& cdata, TID& nextChangePoint)
{
    Init();
    // case one: calibration doesn't exist at all
    if (!dataBase->Has(calibrationID))
        return false;

    // case two: eventID lies inside data, then return it

    return dataBase->GetItem(calibrationID,eventID,cdata,nextChangePoint);
}


const list<TID> DataManager::GetChangePoints(const string& calibrationID)
{
    Init();
    return dataBase->GetChangePoints(calibrationID);
}

uint32_t DataManager::GetNumberOfCalibrations()
{
    Init();
    return dataBase->GetKeys().size();
}

uint32_t DataManager::GetNumberOfDataPoints(const string& calibrationID)
{
    Init();
    return dataBase->GetNumberOfDataPoints(calibrationID);
}
