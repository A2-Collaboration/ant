#include "DataManager.h"

//ant
#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/memory.h"
#include "base/interval.h"
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
    dataBase = std_ext::make_unique<DataBase>(calibrationDataFolder);
}

DataManager::DataManager(const string& calibrationDataFolder_):
    calibrationDataFolder(calibrationDataFolder_)
{}

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
    return dataBase->GetItem(calibrationID,eventID,cdata,nextChangePoint);
}

uint32_t DataManager::GetNumberOfCalibrationIDs()
{
    Init();
    return dataBase->GetKeys().size();
}

uint32_t DataManager::GetNumberOfDataItems(const string& calibrationID)
{
    Init();
    return dataBase->GetNumberOfDataItems(calibrationID);
}
