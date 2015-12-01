#include "DataBase.h"

#include "base/WrapTFile.h"
#include "base/interval.h"
#include "base/Logger.h"

#include "base/std_ext/misc.h" // begins_with
#include "base/std_ext/system.h"

#include "TTree.h"


using namespace std;
using namespace ant;
using namespace ant::calibration;

DataBase::DataBase(const string calibrationDataFolder_):
    calibrationDataFolder(calibrationDataFolder_)
{

}

//bool DataBase::ReadFromFile(const std::string& filename)
//{
//    try {
//        WrapTFileInput dataFile;
//        dataFile.OpenFile(filename);

//        for( TTree* calibtree: dataFile.GetListOf<TTree>())
//        {
//            string tname(calibtree->GetName());
//            if(!std_ext::begins_with(tname, cm_treename_prefix))
//                continue;
//            const TCalibrationData* cdata = nullptr;
//            calibtree->SetBranchAddress(cm_branchname.c_str(),&cdata);
//            for (Long64_t entry = 0; entry < calibtree->GetEntries(); ++entry)
//            {
//                calibtree->GetEntry(entry);
//                dataMap[cdata->CalibrationID].Data.push_back(*cdata); //emplace???
//            }
//        }
//        return true;
//    } catch (...) {
//        LOG(WARNING) << "Cannot open " << filename;
//        return false;
//    }

//}

//bool DataBase::ReadFromFolder(const string& folder)
//{
//    // read all root files in given folder
//    bool flag = false;
//    for(const auto& filename : std_ext::system::lsFiles(folder, ".root")) {
//        ReadFromFile(filename);
//        flag = true;
//    }
//    return flag;
//}

//void DataBase::WriteToTree(WrapTFileOutput& file, const DataMap_t::value_type& calibration) const
//{
//    string tname = cm_treename_prefix + calibration.first;

//    TTree* currentTree = file.CreateInside<TTree>(tname.c_str(),tname.c_str());
//    const TCalibrationData* cdataptr = nullptr;
//    currentTree->Branch(cm_branchname.c_str(), std::addressof(cdataptr));

//    for(const auto& cdata : calibration.second.Data) {
//        cdataptr = std::addressof(cdata);
//        currentTree->Fill();
//    }
//}

//void DataBase::WriteToFile(const std::string& filename) const
//{
//    // write everything into one single file
//    WrapTFileOutput file(filename,
//                         WrapTFileOutput::mode_t::recreate);
//    // loop over map and write a new tree for each calibrationID
//    for(const auto& calibration: dataMap)
//    {
//        WriteToTree(file, calibration);
//    }
//}

//void DataBase::WriteToFolder(const string& folder) const
//{
//    // we cannot use a simple Update mode of the file
//    // since we want to recreate the files from the database
//    // (items might have been deleted...)

//    map<string, shared_ptr<WrapTFileOutput>> files_by_filename;
//    map<string, shared_ptr<WrapTFileOutput>> files_by_id;

//    for(const auto& it_calibration : dataMap)
//    {
//        const string& calibrationID = it_calibration.first;
//        if(calibrationID.empty()) {
//            LOG(WARNING) << "Found empty calibrationID, ignored.";
//            continue;
//        }
//        if(!it_calibration.second.Modified)
//            continue;
//        auto pos = calibrationID.find_first_of('/');
//        string filename = folder + "/" +
//                          calibrationID.substr(0, pos)
//                          + ".root";
//        auto it_file = files_by_filename.find(filename);
//        if(it_file == files_by_filename.end()) {
//            auto file = make_shared<WrapTFileOutput>(filename,
//                                                     WrapTFileOutput::mode_t::recreate);

//            files_by_filename.insert(make_pair(filename, file));
//            files_by_id.insert(make_pair(calibrationID, file));
//        }
//        else {
//            files_by_id.insert(make_pair(calibrationID, it_file->second));
//        }

//    }

//    for(const auto& it_map : files_by_id)
//    {
//        WriteToTree(*it_map.second, *dataMap.find(it_map.first));
//    }
//}

//bool DataBase::Has(const string& calibrationID) const
//{
//    //return dataMap.find(calibrationID) != dataMap.end();
//}



std::list<string> DataBase::GetKeys() const
{
//    list<string> keys;
//    for (const auto& entry: dataMap){
//        keys.emplace_back(entry.first);
//    }
//    return keys;
}

//uint32_t DataBase::getDepth(const TID& tid, const string& calibrationID) const
//{
//    uint32_t current_depth = 0;
//    auto& calibPairs = dataMap.at(calibrationID).Data;

//    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
//    {
//        interval<TID> idint(rit->FirstID,rit->LastID);
//        if (idint.Contains(tid))
//            return current_depth;
//        current_depth++;
//    }
//    return current_depth;
//}

//bool DataBase::isValid(const TID& tid, const string& calibrationID, const uint32_t& depth) const
//{
//    return (depth <= getDepth(tid,calibrationID));
//}



uint32_t DataBase::GetNumberOfDataItems(const string& calibrationID) const
{
//    auto it = dataMap.find(calibrationID);
//    if(it == dataMap.end())
//        return 0;
//    return it->second.Data.size();
}

std::set<DataBase::Range_t> DataBase::getRanges(const string& calibrationID) const
{
    set<Range_t> ranges;
    for(auto rangedir : std_ext::system::lsFiles(calibrationDataFolder+"/"+calibrationID+"/DataRanges")) {
        /// \todo parse folder names
        cout << rangedir << endl;
    }
    return ranges;
}

bool DataBase::loadFile(const string& filename, TCalibrationData& cdata) const
{
    WrapTFileInput dataFile;
    try {
        dataFile.OpenFile(filename);
    }
    catch(...) {
        VLOG(5) << "Cannot open file for reading: " << filename;
        return false;
    }
    return dataFile.GetObjectClone("cdata", cdata);
}

bool DataBase::GetItem(const string& calibrationID,
                       const TID& currentPoint,
                       TCalibrationData& theData,
                       TID& nextChangePoint) const
{
    // try to find it in the ranges
    const auto ranges = getRanges(calibrationID);
    const auto it = ranges.find(currentPoint);
    /// \todo correctly consider ranges, watch out for default data ranges...

    // not found in ranges, so try default data
    if(loadFile(calibrationDataFolder+"/"+calibrationID+"/DataDefault/current", theData)) {
        /// \todo figure out the nextChangePoint
        VLOG(5) << "Loaded default data for " << calibrationID << " for changepoint " << currentPoint;
        return true;
    }

    // nothing found at all
    return false;
}

void DataBase::AddItem(const TCalibrationData& cdata, mode_t mode)
{

}



