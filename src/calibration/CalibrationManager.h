#pragma once

//ant
#include "tree/TCalibrationData.h"
#include "tree/TDataRecord.h"

//ROOT
#include "TFile.h"
#include "TTree.h"

//std
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

namespace ant
{


class CalibrationManager
{
private:
    std::string dataFileName;
    std::vector<TCalibrationData> dataBase;

public:
    CalibrationManager(const std::string& DataFileName):
    dataFileName(DataFileName)
    {
        TFile* dataFile;
        TTree* dataSets;

        dataFile = new TFile(dataFileName.c_str(),"READ");
        if ( dataFile->IsOpen() )
        {
            dataFile->GetObject("dataSets", dataSets);
            if (dataSets != nullptr){
                TCalibrationData* cdata;
                dataSets->SetBranchAddress("dataset",&cdata);
                for ( Long64_t entry = 0; entry < dataSets->GetEntries(); ++entry){
                    dataSets->GetEntry(entry);
                    Add(*cdata);
                }
            }

            dataFile->Close();
            delete(dataFile);
            delete(dataSets);
        }
    }

    void Add(const TCalibrationData& data)
    {
        dataBase.push_back(data);
    }

    const TCalibrationData GetData(const std::string& setupID, const TID& eventID) const
    {
        std::cout << "TODO: add data" << std::endl;
        return TCalibrationData();
    }

    const std::vector<TID> GetChangePoints(const std::string& setupID) const
    {
        return {TID()};
    }

};

} //namespace ant
