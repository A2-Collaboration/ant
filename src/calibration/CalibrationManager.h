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
    TFile* dataFile;
    TTree* dataSets;

public:
    CalibrationManager(const std::string& dataFileName)
    {

        dataFile = new TFile(dataFileName.c_str(),"UPDATE");
        if ( !dataFile->IsOpen() )
            throw false;

        dataSets = dataFile->GetObject("dataSets", dataSets);
        if (dataSets == nullptr)
            dataSets = new TTree("dataSets");

    }

    void Add(const TCalibrationData& data)
    {
        /// TODO: Sanity checks like eg.
        ///  valid setupID?
        ///  only one specific detector per beamtime

    }

    const TCalibrationData& GetData(const std::string& setupID, const TID& eventID) const
    {
        std::cout << "TODO: add data" << endl;
        return TCalibrationData();
    }

    const std::vector<TID> GetChangePoints(const std::string& setupID) const
    {
        return {TID()};
    }

};

} //namespace ant
