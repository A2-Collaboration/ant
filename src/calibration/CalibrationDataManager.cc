#include "CalibrationDataManager.h"

//ant
#include "base/WrapTFile.h"
#include "base/Logger.h"

//ROOT
#include "TTree.h"
#include "TError.h"
#include "TDirectory.h"

//std
#include <algorithm>
#include <iostream>

using namespace std;
using namespace ant;

CalibrationDataManager::CalibrationDataManager(const string& DataFileName):
        cm_treename_prefix("calibration-"),
        cm_branchname("cdata"),
        dataFileName(DataFileName),
        changedDataBase(false)
{
}

void CalibrationDataManager::readDataBase()
{
    try {
        WrapTFile dataFile(dataFileName,
                           WrapTFile::mode_t::read);

        for( TTree* calibtree: dataFile.GetListOf<TTree>())
        {
            const TCalibrationData* cdata = nullptr;
            calibtree->SetBranchAddress(cm_branchname.c_str(),&cdata);
            for (Long64_t entry = 0; entry < calibtree->GetEntries(); ++entry)
            {
                calibtree->GetEntry(entry);
                Add(*cdata);
            }
        }
    } catch (...) {
        VLOG(3) << "Cannot open " << dataFileName << ": Starting with empty database.";
    }

}

void CalibrationDataManager::writeDataBase()
{
    lazyInit();
    WrapTFile file(dataFileName,
                   WrapTFile::mode_t::recreate);
    // loop over map and write a new tree for each calibrationID
    for (auto& calibration: *dataBase)
    {
        string tname = cm_treename_prefix + calibration.first;

        TTree* currentTree = file.CreateInside<TTree>(tname.c_str(),tname.c_str());
        const TCalibrationData* cdataptr = nullptr;
        currentTree->Branch(cm_branchname.c_str(),&cdataptr);

        for( auto& cdata: calibration.second){
            cdataptr = &cdata;
            currentTree->Fill();
        }
    }
}

void CalibrationDataManager::lazyInit()
{
    if ( dataBase == nullptr )
    {
        dataBase = std_ext::make_unique<dataMap>();
        readDataBase();
    }
}

uint32_t CalibrationDataManager::getDepth(const TID& tid, const string& calibrationID) const
{
    uint32_t current_depth = 0;
    auto& calibPairs = dataBase->at(calibrationID);

    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
    {
        interval<TID> idint(rit->FirstID,rit->LastID);
        if (idint.Contains(tid))
            return current_depth;
        current_depth++;
    }
    return current_depth;
}

bool CalibrationDataManager::GetData(const string& calibrationID, const TID& eventID, TCalibrationData& cdata)
{
    lazyInit();
    //case one: calibration doesn't exist
    if ( dataBase->count(calibrationID) == 0)
        return false;

    //case two: calibration exists
    auto& calibPairs = dataBase->at(calibrationID);
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


const list<TID> CalibrationDataManager::GetChangePoints(const string& calibrationID)
{
    lazyInit();
    if ( dataBase->count(calibrationID) == 0)
        return {};

    uint32_t depth = 0;
    list<TID> ids;


    auto& calibPairs = dataBase->at(calibrationID);

    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
    {
        //changepoint is one after the last element;
        auto inclastID(rit->LastID);
        ++inclastID;

        if (isValid(rit->FirstID,calibrationID,depth) )
            ids.push_back(rit->FirstID);
        if (isValid(inclastID,calibrationID,depth) )
            ids.push_back(inclastID);

        depth++;
    }
    ids.sort();
    return ids;
}

uint32_t CalibrationDataManager::GetNumberOfDataPoints(const string& calibrationID)
{
    lazyInit();
    try
    {
        return dataBase->at(calibrationID).size();
    }
    catch (out_of_range)
    {
        return 0;
    }
}

bool CalibrationDataManager::GetIDRange(const string& calibrationID, interval<TID>& IDinterval)
{
    lazyInit();
    if (dataBase->count(calibrationID) == 0)
        return false;

    auto& data = dataBase->at(calibrationID);

    IDinterval.Start() = data.front().FirstID;
    IDinterval.Stop()  = data.front().LastID;

    for (auto& entry: data)
    {
        if (entry.FirstID < IDinterval.Start())
                IDinterval.Start() = entry.FirstID;
        if (entry.LastID  > IDinterval.Stop())
                IDinterval.Stop() = entry.LastID;
    }

    return true;
}

bool CalibrationDataManager::GetLastEntry(const std::string& calibrationID, TCalibrationData& cdata)
{
    lazyInit();
    if (dataBase->count(calibrationID) == 0)
        return false;

    auto& data = dataBase->at(calibrationID);

    cdata = data.back();

    return true;
}
