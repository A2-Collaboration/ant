#pragma once

//ant
#include "tree/TCalibrationData.h"
#include "tree/TDataRecord.h"
#include "base/interval.h"
#include "base/WrapTFile.h"

//
#include "base/Logger.h"

//ROOT
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"
#include "TList.h"

//std
#include <memory>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace ant
{


class CalibrationManager
{
private:
    const std::string cm_treename_prefix;
    const std::string cm_branchname;
    std::string dataFileName;

    std::map<std::string,std::vector<TCalibrationData>> dataBase;

    /**
     * @brief isValid tests if the given id is a valid changepoint ( no newer data exists )
     * @param tid      queried event id
     * @param setupID  calibration id
     * @param depth    depth(distance from last calibration iteration) of given data point
     * @return valid or not
     */
    bool isValid(const TID& tid, const std::string& setupID, const std::uint32_t& depth) const
    {
        return (depth <= getDepth(tid,setupID));
    }

    /**
     * @brief getDepth returns the distance in steps to the latest calibration iteration
     * @param tid      event id
     * @param setupID  calibration id
     * @return depth
     */
    std::uint32_t getDepth(const TID& tid, const std::string& setupID) const
    {
        std::uint32_t current_depth = 0;
        auto& calibPairs = dataBase.at(setupID);

        for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
        {
            interval<TID> idint(rit->FirstID,rit->LastID);
            if (idint.Contains(tid))
                return current_depth;
            current_depth++;
        }
        return current_depth;
    }

    /**
     * @brief finish takes care of rewriting the data to the tree
     */
    void finish() const
    {
        WrapTFile file(dataFileName);
        std::vector<TTree*> treeBuffer;
        // loop over map and write a new tree for each setupID
        for (auto& calibration: dataBase)
        {
            std::string tname = cm_treename_prefix + calibration.first;
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

public:
    CalibrationManager(const std::string& DataFileName):
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
                std::cerr << "no keys  in file " << dataFileName << std::endl;
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
                    std::string treename(calibtree->GetName());
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

    ~CalibrationManager()
    {
        finish();
    }

    void Add(const TCalibrationData& data)
    {
        dataBase[data.SetupID].push_back(data);
    }

    ///
    /// \brief GetData Query the calibration database for specific TID
    /// \param setupID Calibration ID
    /// \param eventID event ID
    /// \param cdata   Reference to a TCalibrationData, data will be writter here
    /// \return true if valid data was found
    ///
    bool GetData(const std::string& setupID, const TID& eventID, TCalibrationData& cdata) const
    {
        //case one: calibration doesn't exist
        if ( dataBase.count(setupID) == 0)
            return false;

        //case two: calibration exists
        auto& calibPairs = dataBase.at(setupID);
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

    const std::vector<TID> GetChangePoints(const std::string& setupID) const
    {
        if ( dataBase.count(setupID) == 0)
            return {};

        std::uint32_t depth = 0;
        std::vector<TID> ids;


        auto& calibPairs = dataBase.at(setupID);

        for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
        {
            //changepoint is one after the last element;
            auto inclastID(rit->LastID);
            ++inclastID;

            if (isValid(rit->FirstID,setupID,depth) )
                ids.push_back(rit->FirstID);
            if (isValid(inclastID,setupID,depth) )
                ids.push_back(inclastID);

            depth++;
        }
        std::sort(ids.begin(),ids.end());
        return ids;
    }

    std::uint32_t GetNumberOfCalibrations() const
    {
        return dataBase.size();
    }

    std::uint32_t GetNumberOfDataPoints(const std::string& setupID) const
    {
        try
        {
            return dataBase.at(setupID).size();
        }
        catch (std::out_of_range)
        {
            return 0;
        }

    }

};

} //namespace ant
