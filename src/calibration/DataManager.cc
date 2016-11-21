#include "DataManager.h"

//ant
#include "DataBase.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/memory.h"
#include "base/interval.h"
#include "tree/TCalibrationData.h"

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

DataManager::~DataManager()
{

}

void DataManager::Add(const TCalibrationData& cdata, Calibration::AddMode_t addMode)
{
    if(override_as_default) {
        addMode = Calibration::AddMode_t::AsDefault;
        LOG(INFO) << "Setting Database Add-Mode to Default";
    }

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

TObject*DataManager::GetTObject(const string& calibrationID, const string& objname, const TID& eventID, TID& nextChangePoint)
{
    Init();
    return dataBase->GetTObject(calibrationID,eventID,objname,nextChangePoint);
}

size_t DataManager::GetNumberOfCalibrationIDs()
{
    Init();
    return dataBase->GetCalibrationIDs().size();
}

size_t DataManager::GetNumberOfCalibrationData(const string& calibrationID)
{
    Init();
    return dataBase->GetNumberOfCalibrationData(calibrationID);
}

bool DataManager::GetOverrideToDefault() const
{
    return override_as_default;
}

void DataManager::SetOverrideToDefault(bool v)
{
    override_as_default = v;
}

string DataManager::GetCalibrationDataFolder() const
{
    return calibrationDataFolder;
}

list<string> ant::calibration::DataManager::GetCalibrationIDs()
{
    Init();
    return dataBase->GetCalibrationIDs();
}
