#include "CalibrationDataManager.h"

//ant
#include "base/WrapTFile.h"
#include "base/interval.h"

//ROOT
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"
#include "TList.h"

//std
#include <algorithm>
#include <iostream>

using namespace std;
using namespace ant;

CalibrationDataManager::CalibrationDataManager(const string& DataFileName):
        cm_treename_prefix("calibration-"),
        cm_branchname("cdata"),
        dataFileName(DataFileName)
{

    TFile dataFile(dataFileName.c_str(),"READ");

    if ( dataFile.IsOpen() )
    {
        TList* keys = dataFile.GetListOfKeys();

        if (!keys)
        {
            cerr << "no keys  in file " << dataFileName << endl;
        }
        else
        {
            TTree* calibtree = nullptr;
            TKey*  key  = nullptr;
            const TCalibrationData* cdata = nullptr;
            TIter nextk(keys);

            while ((key = (TKey*)nextk()))
            {
                calibtree = dynamic_cast<TTree*>(key->ReadObj());
                if ( !calibtree )
                    continue;

                // sanity check: is this the tree you're looking for?
                string treename(calibtree->GetName());
                if (treename.find(cm_treename_prefix) != 0)
                    continue;

                calibtree->SetBranchAddress(cm_branchname.c_str(),&cdata);
                for (Long64_t entry = 0; entry < calibtree->GetEntries(); ++entry)
                {
                    calibtree->GetEntry(entry);
                    Add(*cdata);
                }
            }
        }

        dataFile.Close();
    }
}

void CalibrationDataManager::finish() const
{
    WrapTFile file(dataFileName);
    vector<TTree*> treeBuffer;
    // loop over map and write a new tree for each calibrationID
    for (auto& calibration: dataBase)
    {
        string tname = cm_treename_prefix + calibration.first;
        TTree* currentTree = new TTree(tname.c_str(),tname.c_str());
        const TCalibrationData* cdataptr = nullptr;
        currentTree->Branch(cm_branchname.c_str(),&cdataptr);
        for( auto& cdata: calibration.second){
            cdataptr = &cdata;
            currentTree->Fill();
        }
        treeBuffer.push_back(currentTree);
    }
}

uint32_t CalibrationDataManager::getDepth(const TID& tid, const string& calibrationID) const
{
    uint32_t current_depth = 0;
    auto& calibPairs = dataBase.at(calibrationID);

    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
    {
        interval<TID> idint(rit->FirstID,rit->LastID);
        if (idint.Contains(tid))
            return current_depth;
        current_depth++;
    }
    return current_depth;
}

bool CalibrationDataManager::GetData(const string& calibrationID, const TID& eventID, TCalibrationData& cdata) const
{
    //case one: calibration doesn't exist
    if ( dataBase.count(calibrationID) == 0)
        return false;

    //case two: calibration exists
    auto& calibPairs = dataBase.at(calibrationID);
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


const list<TID> CalibrationDataManager::GetChangePoints(const string& calibrationID) const
{
    if ( dataBase.count(calibrationID) == 0)
        return {};

    uint32_t depth = 0;
    list<TID> ids;


    auto& calibPairs = dataBase.at(calibrationID);

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

uint32_t CalibrationDataManager::GetNumberOfDataPoints(const string& calibrationID) const
{
    try
    {
        return dataBase.at(calibrationID).size();
    }
    catch (out_of_range)
    {
        return 0;
    }

}
