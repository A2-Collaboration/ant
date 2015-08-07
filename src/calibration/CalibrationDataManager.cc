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
using namespace ant::calibration;

void DataBase::ReadData(const std::string& filename)
{
    try {
        WrapTFile dataFile(filename,
                           WrapTFile::mode_t::read);

        for( TTree* calibtree: dataFile.GetListOf<TTree>())
        {
            const TCalibrationData* cdata = nullptr;
            calibtree->SetBranchAddress(cm_branchname.c_str(),&cdata);
            for (Long64_t entry = 0; entry < calibtree->GetEntries(); ++entry)
            {
                calibtree->GetEntry(entry);
                DataMap[cdata->CalibrationID].push_back(*cdata); //emplace???
            }
        }
    } catch (...) {
        VLOG(3) << "Cannot open " << filename;
    }

}

void DataBase::WriteData(const std::string& filename) const
{
    WrapTFile file(filename,
                   WrapTFile::mode_t::recreate);
    // loop over map and write a new tree for each calibrationID
    for (auto& calibration: DataMap)
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

void DataManager::lazyInit()
{
    if ( dataBase == nullptr )
    {
        dataBase = std_ext::make_unique<DataBase>();
        dataBase->ReadData(dataFileName);
    }
}

uint32_t DataManager::getDepth(const TID& tid, const string& calibrationID) const
{
    uint32_t current_depth = 0;
    auto& calibPairs = dataBase->DataMap.at(calibrationID);

    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
    {
        interval<TID> idint(rit->FirstID,rit->LastID);
        if (idint.Contains(tid))
            return current_depth;
        current_depth++;
    }
    return current_depth;
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
    if ( dataBase->DataMap.count(calibrationID) == 0)
        return {};

    uint32_t depth = 0;
    list<TID> ids;

    auto& calibPairs = dataBase->DataMap.at(calibrationID);

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

uint32_t DataManager::GetNumberOfDataPoints(const string& calibrationID)
{
    lazyInit();
    try
    {
        return dataBase->DataMap.at(calibrationID).size();
    }
    catch (out_of_range)
    {
        return 0;
    }
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
