#include "DataBase.h"

#include "base/WrapTFile.h"
#include "base/interval.h"


#include "TTree.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;



bool DataBase::ReadData(const std::string& filename)
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
        return true;
    } catch (...) {
        //VLOG(3) << "Cannot open " << filename;
        return false;
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

uint32_t DataBase::getDepth(const TID& tid, const string& calibrationID) const
{
    uint32_t current_depth = 0;
    auto& calibPairs = DataMap.at(calibrationID);

    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
    {
        interval<TID> idint(rit->FirstID,rit->LastID);
        if (idint.Contains(tid))
            return current_depth;
        current_depth++;
    }
    return current_depth;
}

uint32_t DataBase::GetNumberOfDataPoints(const string& calibrationID) const
{
    try
    {
        return DataMap.at(calibrationID).size();
    }
    catch (out_of_range)
    {
        return 0;
    }
}

const std::list<TID> DataBase::GetChangePoints(const string& calibrationID) const
{
    if ( DataMap.count(calibrationID) == 0)
        return {};

    uint32_t depth = 0;
    list<TID> ids;

    auto& calibPairs = DataMap.at(calibrationID);

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
